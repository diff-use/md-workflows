# md-workflows

A Python CLI and Docker image for running the full molecular-dynamics workflow pipeline
(GROMACS + AmberTools + DIALS/cctbx + lunus).

## Docker image layout

The image is built in stages so the architecture-specific GROMACS compile is isolated from
the portable scientific stack — forks targeting other hardware can swap the GROMACS stage
without rebuilding the base:

| Stage | Dockerfile | Contents |
|-------|-----------|----------|
| base    | `Dockerfile.base`    | CUDA 12.6 devel toolchain, micromamba/conda `lunus` env, lunus, ChimeraX. Architecture-neutral. |
| gromacs | `Dockerfile.gromacs` | GROMACS (CUDA, tuned for H100 / AVX-512 by default) + the `md_workflows` package. The consumable image. |

(An Astera-specific `Dockerfile.actl` overlay adds workspace conventions on top of `gromacs`. It
is tracked here so it stays under CI lint coverage, but is built and published only from the
`astera` deployment branch.)

## 1) Build the images

From the project root:

```bash
# 1. Base — architecture-neutral foundation
docker build -f Dockerfile.base -t md-base:local .

# 2. GROMACS + md-workflows (the runnable image).
#    Override the GMX_* build args for non-H100 / non-AVX-512 hardware, e.g.
#    --build-arg GMX_CUDA_TARGET_SM=80 --build-arg GMX_SIMD=AVX2_256
docker build -f Dockerfile.gromacs --build-arg BASE_IMAGE=md-base:local -t md-gromacs:local .
```

CI builds these stages and pushes versioned tags (derived from the `version` in
`pyproject.toml`) to the Astera Harbor registry; see `.github/workflows/build-images.yml` on
the `astera` branch.

## 2) Start a container

Run interactively, mounting the project directory so inputs/outputs are available on the host:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  --gpus all \
  --name md_container \
  -e HOME=/workspace \
  -v "$(pwd):/workspace" \
  -w /workspace \
  md-gromacs:local \
  bash
```

This registers the CLI entry points from `pyproject.toml`, including `md_workflows.mdmx`.

## 3) Run the full workflow command

Inside the container shell:

```bash
md_workflows.mdmx \
  --param-pdb-id 6B8X \
  --ix 1 \
  --ntomp 26 \
  --resolv-ntmpi 8 \
  --resolv-ntomp 1
```

To see all available flags:

```bash
md_workflows.mdmx --help
```

## Development

Lint and format with [ruff](https://docs.astral.sh/ruff/) (config in `pyproject.toml`):

```bash
pip install '.[dev]'
ruff check md_workflows
ruff format --check md_workflows
```
