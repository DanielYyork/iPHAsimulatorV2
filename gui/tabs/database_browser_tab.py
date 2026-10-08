"""Browse and download files from the structure database."""

import mimetypes
from pathlib import Path

import streamlit as st

from gui.config import STRUCTURE_DATABASE


PREVIEW_BYTES = 64 * 1024
LARGE_DOWNLOAD_BYTES = 200 * 1024 * 1024


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


def render_database_browser_tab(root: Path = STRUCTURE_DATABASE) -> None:
    """Render a folder browser with one-file-at-a-time downloads."""
    st.markdown("## 📁 Structure Database")
    st.write("Browse generated files and download a selected file to your computer.")

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

    shown_folder = Path("structure_database") / relative
    st.code(str(shown_folder))

    if relative != Path(".") and st.button("⬆️ Up one folder", key="database_up"):
        st.session_state[state_key] = str(relative.parent)
        st.rerun()

    if not entries:
        st.info("This folder is empty.")
        return

    filter_text = st.text_input(
        "Filter this folder", key=f"database_filter_{relative}"
    ).casefold()
    visible = [entry for entry in entries if filter_text in entry.name.casefold()]
    if not visible:
        st.info("No files or folders match this filter.")
        return

    names = [entry.name for entry in visible]
    chosen = st.selectbox(
        "Select a folder or file",
        names,
        format_func=lambda name: (
            f"📁 {name}" if (root / relative / name).is_dir() else f"📄 {name}"
        ),
        key=f"database_entry_{relative}",
    )
    selected_relative = relative / chosen

    try:
        selected = database_path(root, selected_relative)
    except (OSError, ValueError) as error:
        st.warning(str(error))
        return

    if selected.is_dir():
        if st.button("Open folder", key="database_open"):
            st.session_state[state_key] = str(selected_relative)
            st.rerun()
        return

    size = selected.stat().st_size
    st.write(f"**File:** `{selected_relative}`")
    st.write(f"**Size:** {_size_label(size)}")

    if size > LARGE_DOWNLOAD_BYTES:
        st.warning(
            "This is a large file. Downloading it through Streamlit may use "
            "substantial server memory."
        )

    mime = mimetypes.guess_type(selected.name)[0] or "application/octet-stream"
    st.download_button(
        "⬇️ Download selected file",
        data=lambda: download_file(root, selected_relative),
        file_name=selected.name,
        mime=mime,
        key=f"database_download_{selected_relative}",
    )

    if size <= PREVIEW_BYTES:
        with selected.open("rb") as file:
            sample = file.read(PREVIEW_BYTES)
        if b"\0" not in sample:
            try:
                preview = sample.decode("utf-8")
            except UnicodeDecodeError:
                pass
            else:
                with st.expander("Preview file"):
                    st.code(preview)
