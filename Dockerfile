# syntax=docker/dockerfile:1
# The explicit package lock is for linux/amd64. Render uses this architecture.
FROM mambaorg/micromamba:2.9.0

COPY --chown=$MAMBA_USER:$MAMBA_USER deploy/conda-linux-64.lock /tmp/conda-linux-64.lock
# @EXPLICIT installs exact archives with checksums: no dependency solve here.
# Keep the package cache: micromamba clean --all fails on AmberTools 26's
# cached bin/amber.conda entry after installation has already completed.
# Upstream: https://github.com/conda-forge/ambertools-feedstock/issues/193
RUN micromamba install --yes --name base --file /tmp/conda-linux-64.lock

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/app:/app/src \
    IPHASIMULATOR_PYTHON=/opt/conda/bin/python \
    AMBERHOME=/opt/conda \
    IPHA_OPENMM_PLATFORM=CPU \
    OPENMM_CPU_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/ipha-matplotlib \
    NUMBA_CACHE_DIR=/tmp/ipha-numba

COPY --chown=$MAMBA_USER:$MAMBA_USER deploy/pip-requirements.txt deploy/pip-constraints.txt /tmp/
RUN python -m pip install --no-cache-dir -c /tmp/pip-constraints.txt -r /tmp/pip-requirements.txt && \
    python -m pip check

USER root
RUN mkdir -p /app /opt/ipha-seed && \
    chown "$MAMBA_USER:$MAMBA_USER" /app /opt/ipha-seed
USER $MAMBA_USER
WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml pha_gui.py streamlit_cloud.py ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ ./src/
COPY --chown=$MAMBA_USER:$MAMBA_USER gui/ ./gui/
COPY --chown=$MAMBA_USER:$MAMBA_USER cluster/ ./cluster/
COPY --chown=$MAMBA_USER:$MAMBA_USER cloud/streamlit_cloud.py cloud/verify.py ./cloud/
COPY --chown=$MAMBA_USER:$MAMBA_USER deploy/ ./deploy/
# Only the chemistry library is shipped; research systems/trajectories are not.
COPY --chown=$MAMBA_USER:$MAMBA_USER structure_database/residue_codes.csv /opt/ipha-seed/residue_codes.csv
COPY --chown=$MAMBA_USER:$MAMBA_USER structure_database/PHA_types/ /opt/ipha-seed/PHA_types/
RUN python -m pip install --no-cache-dir --no-deps --no-build-isolation . && \
    python -m pip check && \
    python deploy/start.py --prepare-only && \
    python cloud/verify.py --output /tmp/ipha-build-check && \
    cp /tmp/ipha-build-check/report.json /app/deploy/build-report.json

# A fresh persistent disk can be root-owned. The launcher sets its top-level
# ownership, then drops to the image's unprivileged user before starting the GUI.
USER root
EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
    CMD ["/opt/conda/bin/python", "/app/deploy/healthcheck.py"]
# Retain micromamba's inherited ENTRYPOINT, which activates the environment.
CMD ["python", "/app/deploy/start.py"]
