"""Browse and download files from the structure database."""

import mimetypes
from datetime import datetime
from pathlib import Path

import streamlit as st

from gui.config import STRUCTURE_DATABASE


PREVIEW_BYTES = 64 * 1024
IMAGE_PREVIEW_BYTES = 10 * 1024 * 1024
LARGE_DOWNLOAD_BYTES = 200 * 1024 * 1024
PAGE_SIZE = 40


def database_path(root: Path, relative: str | Path) -> Path:
    """Resolve a database entry without allowing traversal or symlink escapes."""
    root = Path(root).resolve(strict=True)
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("The selected path is outside the structure database.")

    candidate = root
    for part in relative.parts:
        candidate = candidate / part
        if candidate.is_symlink():
            raise ValueError("Linked files and folders cannot be downloaded.")

    resolved = candidate.resolve(strict=True)
    if not resolved.is_relative_to(root):
        raise ValueError("The selected path is outside the structure database.")
    return resolved


def database_entries(root: Path, relative: str | Path) -> list[Path]:
    """List just the selected folder so large databases remain browsable."""
    folder = database_path(root, relative)
    if not folder.is_dir():
        raise NotADirectoryError(folder)
    return sorted(
        folder.iterdir(),
        key=lambda path: (not path.is_dir(), path.name.casefold()),
    )


def download_file(root: Path, relative: str | Path) -> bytes:
    """Recheck the selected file when a user actually requests its contents."""
    file_path = database_path(root, relative)
    if not file_path.is_file():
        raise ValueError("The selected item is not a file.")
    return file_path.read_bytes()


def _size_label(size: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{size} B"
        size /= 1024
    raise AssertionError("Unreachable")


def _open_folder(relative: Path) -> None:
    st.session_state.structure_database_browser_folder = str(relative)
    st.session_state.structure_database_browser_selected = None
    st.rerun()


def _render_folder_tree(root: Path, relative: Path, entries: list[Path]) -> None:
    """Show root shortcuts, the current path, and nearby folders."""
    st.markdown("#### Folders")
    if st.button("📁 Structure Database", key="database_tree_root", width="stretch"):
        _open_folder(Path("."))

    ancestors = list(reversed(relative.parents)) + [relative]
    for index, ancestor in enumerate(ancestors):
        if ancestor == Path("."):
            continue
        if st.button(
            f"↳ {ancestor.name}",
            key=f"database_tree_ancestor_{index}",
            width="stretch",
            type="primary" if ancestor == relative else "secondary",
        ):
            _open_folder(ancestor)

    st.divider()
    if relative != Path(".") and st.button(
        "⬆️ Parent folder", key="database_up", width="stretch"
    ):
        _open_folder(relative.parent)

    root_folders = [
        path for path in database_entries(root, ".")
        if path.is_dir() and not path.is_symlink()
    ]
    with st.container(height=450):
        st.caption("Top-level folders")
        for path in root_folders:
            if st.button(
                f"📁 {path.name}", key=f"database_root_{path.name}", width="stretch"
            ):
                _open_folder(Path(path.name))

        if relative != Path("."):
            child_folders = [
                path for path in entries if path.is_dir() and not path.is_symlink()
            ]
            if child_folders:
                st.caption("Folders here")
            for path in child_folders[:20]:
                if st.button(
                    f"📁 {path.name}",
                    key=f"database_child_{relative / path.name}",
                    width="stretch",
                ):
                    _open_folder(relative / path.name)
            if len(child_folders) > 20:
                st.caption("More folders are in the file list →")


def _sorted_entries(entries: list[Path], sort_by: str, descending: bool) -> list[Path]:
    """Put directories first, then sort each group by the selected column."""
    def value(path: Path):
        if sort_by == "Size":
            return path.lstat().st_size if path.is_file() or path.is_symlink() else 0
        if sort_by == "Modified":
            return path.lstat().st_mtime
        if sort_by == "Type":
            return path.suffix.casefold(), path.name.casefold()
        return path.name.casefold()

    folders = [path for path in entries if path.is_dir() and not path.is_symlink()]
    folder_set = set(folders)
    files = [path for path in entries if path not in folder_set]
    return sorted(folders, key=value, reverse=descending) + sorted(
        files, key=value, reverse=descending
    )


def _render_file_list(root: Path, relative: Path, entries: list[Path]) -> None:
    st.markdown("#### Files and folders")
    parts = ["structure_database", *relative.parts] if relative != Path(".") else ["structure_database"]
    with st.container(horizontal=True):
        for index, name in enumerate(parts):
            destination = Path(*parts[1:index + 1]) if index else Path(".")
            if st.button(
                name, key=f"database_crumb_{relative}_{index}"
            ):
                _open_folder(destination)

    filter_column, sort_column, order_column = st.columns([3, 1.3, 1.2])
    with filter_column:
        query = st.text_input(
            "Search this folder", key=f"database_filter_{relative}"
        ).casefold()
    with sort_column:
        sort_by = st.selectbox(
            "Sort by", ["Name", "Type", "Size", "Modified"],
            key=f"database_sort_{relative}",
        )
    with order_column:
        descending = st.toggle("Descending", key=f"database_reverse_{relative}")

    visible = _sorted_entries(
        [entry for entry in entries if query in entry.name.casefold()],
        sort_by,
        descending,
    )
    st.caption(f"{len(visible)} of {len(entries)} items")
    if not visible:
        st.info("No files or folders match this search." if entries else "This folder is empty.")
        return

    total_pages = (len(visible) + PAGE_SIZE - 1) // PAGE_SIZE
    page_key = f"database_page_{relative}"
    page = min(st.session_state.get(page_key, 0), total_pages - 1)
    st.session_state[page_key] = page
    if total_pages > 1:
        previous, page_label, next_page = st.columns([1, 2, 1])
        with previous:
            if st.button("← Previous", disabled=page == 0, key=f"database_previous_{relative}"):
                st.session_state[page_key] = page - 1
                st.rerun()
        with page_label:
            st.caption(f"Page {page + 1} of {total_pages}")
        with next_page:
            if st.button("Next →", disabled=page == total_pages - 1, key=f"database_next_{relative}"):
                st.session_state[page_key] = page + 1
                st.rerun()

    header = st.columns([5, 1.4, 1.2, 1.8])
    for column, label in zip(header, ("Name", "Type", "Size", "Modified")):
        column.markdown(f"**{label}**")
    st.divider()

    selected_name = st.session_state.get("structure_database_browser_selected")
    with st.container(height=450):
        for path in visible[page * PAGE_SIZE:(page + 1) * PAGE_SIZE]:
            is_link = path.is_symlink()
            is_folder = path.is_dir() and not is_link
            row = st.columns([5, 1.4, 1.2, 1.8], vertical_alignment="center")
            relative_path = relative / path.name
            with row[0]:
                if st.button(
                    f"{'📁' if is_folder else '🔗' if is_link else '📄'} {path.name}",
                    key=f"database_row_{relative_path}",
                    width="stretch",
                    type="primary" if str(relative_path) == selected_name else "secondary",
                    disabled=is_link,
                    help="Open folder" if is_folder else "Select file to preview and download",
                ):
                    if is_folder:
                        _open_folder(relative_path)
                    st.session_state.structure_database_browser_selected = str(relative_path)
                    st.rerun()
            row[1].write("Folder" if is_folder else "Link" if is_link else (path.suffix or "File"))
            try:
                details = path.stat() if not is_link else path.lstat()
            except OSError:
                row[2].write("—")
                row[3].write("—")
            else:
                row[2].write("—" if is_folder else _size_label(details.st_size))
                row[3].write(datetime.fromtimestamp(details.st_mtime).strftime("%Y-%m-%d"))


def _render_selection(root: Path, relative: Path) -> None:
    selected_name = st.session_state.get("structure_database_browser_selected")
    if not selected_name:
        st.info("Select a file to preview it and download a copy.")
        return

    selected_relative = Path(selected_name)
    if selected_relative.parent != relative:
        st.session_state.structure_database_browser_selected = None
        return
    try:
        selected = database_path(root, selected_relative)
        if not selected.is_file():
            raise ValueError("The selected item is not a file.")
        size = selected.stat().st_size
    except (OSError, ValueError) as error:
        st.warning(f"The selected file is unavailable: {error}")
        return

    st.markdown(f"#### 📄 {selected.name}")
    st.caption(f"{selected_relative} · {_size_label(size)}")
    if size > LARGE_DOWNLOAD_BYTES:
        st.warning(
            "This is a large file. Downloading it through Streamlit may use "
            "substantial server memory."
        )

    mime = mimetypes.guess_type(selected.name)[0] or "application/octet-stream"
    st.download_button(
        "⬇️ Download file",
        data=lambda: download_file(root, selected_relative),
        file_name=selected.name,
        mime=mime,
        key=f"database_download_{selected_relative}",
        type="primary",
    )

    if selected.suffix.casefold() in {".png", ".jpg", ".jpeg", ".webp"} and size <= IMAGE_PREVIEW_BYTES:
        with st.expander("Preview image", expanded=True):
            st.image(selected.read_bytes(), width="stretch")
        return

    if size <= PREVIEW_BYTES:
        with selected.open("rb") as file:
            sample = file.read(PREVIEW_BYTES)
        if b"\0" not in sample:
            try:
                preview = sample.decode("utf-8")
            except UnicodeDecodeError:
                pass
            else:
                with st.expander("Preview file", expanded=True):
                    st.code(preview)


def render_database_browser_tab(root: Path = STRUCTURE_DATABASE) -> None:
    """Render a two-pane browser with one-file-at-a-time downloads."""
    st.markdown("## 📁 Structure Database")
    st.write("Browse folders, inspect files, and download a copy to your computer.")

    root = Path(root)
    if not root.is_dir():
        st.warning("The structure database is not available.")
        return

    state_key = "structure_database_browser_folder"
    relative = Path(st.session_state.get(state_key, "."))
    try:
        entries = database_entries(root, relative)
    except (OSError, ValueError):
        relative = Path(".")
        st.session_state[state_key] = "."
        entries = database_entries(root, relative)

    left, right = st.columns([1, 3.2], gap="medium")
    with left, st.container(border=True):
        _render_folder_tree(root, relative, entries)
    with right, st.container(border=True):
        _render_file_list(root, relative, entries)
        st.divider()
        _render_selection(root, relative)
