# ==================== CUDA 12.6 (single-stage, devel) ====================
# CUDA 12.6 *devel* image: nvcc + full CUDA Toolkit, kept for the final runtime so GROMACS builds
# with -DGMX_GPU=CUDA AND future CUDA-dependent tooling can compile inside the container.
# Host needs the NVIDIA driver + NVIDIA Container Toolkit at run time (e.g. docker run --gpus all).
#   docker build --build-arg CUDA_IMAGE_TAG=12.6.3 -t diffuseproject/md:gpu .
ARG CUDA_IMAGE_TAG=12.6.3
FROM nvidia/cuda:${CUDA_IMAGE_TAG}-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

# ---------- system packages (build tools + runtime deps) ----------
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        wget \
        git \
        bzip2 \
        coreutils \
        rsync \
        bc \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# ---------- install micromamba ----------
# Install envs under /opt so shebangs (#!/opt/micromamba/...) match the runtime layout.
RUN curl -L micro.mamba.pm/install.sh -o /tmp/micromamba_install.sh \
    && printf '\n\n\n\n' | bash /tmp/micromamba_install.sh \
    && rm /tmp/micromamba_install.sh \
    && mkdir -p /opt/micromamba/bin \
    && cp /root/.local/bin/micromamba /opt/micromamba/bin/micromamba

ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV MAMBA_EXE=/opt/micromamba/bin/micromamba
ENV PATH="/opt/micromamba/bin:${PATH}"

# ---------- conda environment (inline of lunus.yaml) ----------
RUN cat > /tmp/lunus.yaml <<'YAML'
name: lunus
channels:
  - conda-forge
dependencies:
  - python =3.10
  - ambertools =24.8
  - dials
  - cctbx-base
  - nexusformat =2.0.2
  - pandas =2.2.3
  - scipy =1.14.1
  - numexpr =2.14.1
  - joblib =1.5.3
  - matplotlib =3.10.0
  - ipykernel =7.2.0
  - ipython =8.32.0
  - jupyterlab =4.5.4
  - nb_conda_kernels =2.5.1
  - gtk3 =3.24.43
  - xarray =2025.01.0
  - nexpy =2.0.1
  - scons =4.10.1
  - git
  - vim
  - curl
  - scp =0.15.0
  - mdtraj =1.10.3
  - openmpi <5
  - mpi4py
  - openssh
  - cmake =3.31.2
  - awscli
  - gnuplot =5.4.10
YAML

RUN $MAMBA_EXE create -y -f /tmp/lunus.yaml && rm /tmp/lunus.yaml

ARG MAMBA_ENV=/opt/micromamba/envs/lunus
ENV PATH="${MAMBA_ENV}/bin:${PATH}"
ENV CONDA_PREFIX="${MAMBA_ENV}"

# ---------- pip packages ----------
RUN pip install --no-cache-dir git+https://github.com/ando-lab/mdx2.git

# ---------- GROMACS (CUDA build, targeting H100 / sm_90) ----------
# GMX_SIMD is pinned to AVX_512 (Voltage Park Xeon Platinum supports it) instead of letting CMake
# auto-detect from the build host: under QEMU emulation detection falls back to SSE4.1, which would
# cripple CPU-side kernels. Pinning makes the CPU SIMD deployment-correct regardless of build host.
RUN set -ex \
    && d=$(mktemp -d) \
    && cd "$d" \
    && wget https://ftp.gromacs.org/gromacs/gromacs-2025.2.tar.gz \
    && tar xfz gromacs-2025.2.tar.gz \
    && cd gromacs-2025.2 \
    && mkdir build && cd build \
    && cmake .. \
        -DGMX_BUILD_OWN_FFTW=ON \
        -DGMX_GPU=CUDA \
        -DCUDAToolkit_ROOT=/usr/local/cuda \
        -DGMX_CUDA_TARGET_SM=90 \
        -DGMX_SIMD=AVX_512 \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf "$d"

# ---------- lunus ----------
RUN mkdir -p /opt/packages \
    && cd /opt/packages \
    && git clone https://github.com/lanl/lunus \
    && cd lunus \
    && scons enable-openmp=True

# ---------- cleanup: remove build-only packages + caches ----------
RUN $MAMBA_EXE remove -n lunus -y scons cmake \
    && $MAMBA_EXE clean -afy \
    && find /opt/micromamba -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null; \
    find /opt/micromamba -name "*.pyc" -delete 2>/dev/null; \
    rm -rf /opt/packages/lunus/.git; \
    true

# ---------- ChimeraX ----------
ARG CHIMERAX_URL="https://www.cgl.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py?file=current/ubuntu-22.04/chimerax-daily.deb"
RUN apt-get update \
    && curl -s -c /tmp/cx_cookies -d "choice=Accept" "${CHIMERAX_URL}" \
       | grep -oP 'url=\K[^"]*' > /tmp/cx_redirect \
    && curl -s -b /tmp/cx_cookies -o /tmp/chimerax.deb \
       "https://www.cgl.ucsf.edu$(cat /tmp/cx_redirect)" \
    && apt-get install -y /tmp/chimerax.deb \
    && rm -f /tmp/chimerax.deb /tmp/cx_cookies /tmp/cx_redirect \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /home/dev/.config/ChimeraX

# ---------- shell init (minimal) ----------
# Runs as root by default (see end of file); downstream images (e.g. Dockerfile.astera) own the
# general user environment. Put the micromamba hook in the *global* bashrc so `micromamba activate`
# works for root and for any UID supplied via `docker run --user ...` (skel only covers new users).
RUN printf '\n# Enable `micromamba activate` in interactive shells\neval "$(micromamba shell hook --shell bash)"\n' >> /etc/bash.bashrc

# ---------- md-workflows ----------
# Ship md-workflows in the lunus env so Hub users need not pip install / extend PATH.
# Placed late so code edits don't invalidate the expensive conda/GROMACS layers.
COPY pyproject.toml /opt/md-workflows/pyproject.toml
COPY md_workflows /opt/md-workflows/md_workflows
RUN /opt/micromamba/envs/lunus/bin/python -m pip install --no-cache-dir /opt/md-workflows

ENV PATH="/opt/micromamba/bin:/opt/micromamba/envs/lunus/bin:/usr/local/gromacs/bin:${PATH}"

# ---------- GPU runtime metadata (placed late so it doesn't bust the build cache) ----------
# Do NOT auto-claim GPUs. The base nvidia/cuda image sets NVIDIA_VISIBLE_DEVICES=all, which under
# Kubernetes (NVIDIA device plugin / GPU Operator) overrides per-pod GPU isolation and exposes every
# node GPU regardless of resource requests. Override to "void" so GPUs are granted only at run time:
# `docker run --gpus ...` and the K8s device plugin both set NVIDIA_VISIBLE_DEVICES themselves.
ENV NVIDIA_VISIBLE_DEVICES=void
# Capabilities to mount when a GPU *is* granted (harmless when none is).
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility

# Run as root by default. No baked non-root user: downstream images (Dockerfile.astera) run as root,
# and the standalone README workflow overrides identity with `docker run --user "$(id -u):$(id -g)"`.
# /opt artifacts stay root-owned at default perms, so they remain readable/executable by any UID.
# A dedicated user can be added later if a use case needs one.
WORKDIR /workspace
SHELL ["/bin/bash", "-c"]
CMD ["bash"]
