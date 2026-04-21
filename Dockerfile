# ==================== BUILDER ====================
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        wget \
        git \
        bzip2 \
    && rm -rf /var/lib/apt/lists/*

# ---------- install micromamba ----------
RUN curl -L micro.mamba.pm/install.sh -o /tmp/micromamba_install.sh \
    && printf '\n\n\n\n' | bash /tmp/micromamba_install.sh \
    && rm /tmp/micromamba_install.sh

ENV MAMBA_ROOT_PREFIX=/root/micromamba
ENV MAMBA_EXE=/root/.local/bin/micromamba
ENV PATH="/root/.local/bin:${PATH}"

# ---------- conda environment (inline of lunus.yaml) ----------
RUN cat > /tmp/lunus.yaml <<'YAML'
name: lunus
channels:
  - conda-forge
dependencies:
  - python =3.10
  - ambertools =24.8
  - dials
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

ARG MAMBA_ENV=/root/micromamba/envs/lunus
ENV PATH="${MAMBA_ENV}/bin:${PATH}"
ENV CONDA_PREFIX="${MAMBA_ENV}"

# ---------- pip packages ----------
RUN pip install --no-cache-dir git+https://github.com/ando-lab/mdx2.git

# ---------- GROMACS (inline of install_gromacs.sh) ----------
RUN set -ex \
    && d=$(mktemp -d) \
    && cd "$d" \
    && wget https://ftp.gromacs.org/gromacs/gromacs-2025.2.tar.gz \
    && tar xfz gromacs-2025.2.tar.gz \
    && cd gromacs-2025.2 \
    && mkdir build && cd build \
    && cmake .. \
        -DGMX_BUILD_OWN_FFTW=ON \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf "$d"

# ---------- lunus ----------
RUN mkdir -p /root/packages \
    && cd /root/packages \
    && git clone https://github.com/lanl/lunus \
    && cd lunus \
    && scons enable-openmp=True

# ---------- cleanup: remove build-only packages + caches ----------
RUN $MAMBA_EXE remove -n lunus -y scons cmake \
    && $MAMBA_EXE clean -afy \
    && find /root/micromamba -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null; \
    find /root/micromamba -name "*.pyc" -delete 2>/dev/null; \
    rm -rf /root/packages/lunus/.git; \
    true

# ==================== FINAL ====================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

ARG CHIMERAX_URL="https://www.cgl.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py?file=current/ubuntu-22.04/chimerax-daily.deb"
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        coreutils \
        rsync \
        bc \
        bzip2 \
        curl \
        libgomp1 \
    && curl -s -c /tmp/cx_cookies -d "choice=Accept" "${CHIMERAX_URL}" \
       | grep -oP 'url=\K[^"]*' > /tmp/cx_redirect \
    && curl -s -b /tmp/cx_cookies -o /tmp/chimerax.deb \
       "https://www.cgl.ucsf.edu$(cat /tmp/cx_redirect)" \
    && apt-get install -y /tmp/chimerax.deb \
    && rm -f /tmp/chimerax.deb /tmp/cx_cookies /tmp/cx_redirect \
    && apt-get purge -y curl \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# ---------- bashrc (inline of bashrc_new) ----------
RUN cat > /root/.bashrc <<'BASHRC'
# If not running interactively, don't do anything
#case $- in
#    *i*) ;;
#      *) return;;
#esac

HISTCONTROL=ignoreboth
shopt -s histappend
HISTSIZE=1000
HISTFILESIZE=2000
shopt -s checkwinsize

[ -x /usr/bin/lesspipe ] && eval "$(SHELL=/bin/sh lesspipe)"

if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
        color_prompt=yes
    else
        color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
fi
unset color_prompt force_color_prompt

case "$TERM" in
xterm*|rxvt*)
    PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
    ;;
*)
    ;;
esac

if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi

eval "$(micromamba shell hook --shell bash)"
BASHRC

# ---------- copy artifacts from builder ----------
COPY --from=builder /root/.local/bin/micromamba /root/.local/bin/micromamba
COPY --from=builder /root/micromamba /root/micromamba
COPY --from=builder /usr/local/gromacs /usr/local/gromacs
COPY --from=builder /root/packages /root/packages

ENV MAMBA_ROOT_PREFIX=/root/micromamba
ENV MAMBA_EXE=/root/.local/bin/micromamba
ENV PATH="/root/.local/bin:/root/micromamba/envs/lunus/bin:/usr/local/gromacs/bin:${PATH}"
ENV CONDA_PREFIX=/root/micromamba/envs/lunus

WORKDIR /root
SHELL ["/bin/bash", "-c"]
CMD ["bash"]
