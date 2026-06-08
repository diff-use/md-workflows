# md-workflows

This project includes a `Dockerfile` and a Python CLI entry point for running the full MD workflow pipeline.

## Running on ACTL

The Astera overlay image is published as:

```bash
harbor.astera.sh/library/md-workflows:0.0.2-actl-2026-06-08
```

It is available in ACTL as the `md-workflows` image alias for the `diffuse` namespace. From this checkout:

```bash
actl pod profiles -n diffuse
actl pod up md-workflows --profile single --image md-workflows --pvc-size 100Gi -n diffuse --yes
```

The selected diffuse profile auto-mounts the shared volume at `/mnt/diffuse-shared`. Keep large inputs, trajectories, and outputs there rather than in the synced source checkout. ACTL syncs this repository to `/home/dev/workspace`; the image also puts `/home/dev/workspace` on `PYTHONPATH`, so local source edits override the baked `/opt/md-workflows` package.

Inside the ACTL shell:

```bash
cd /mnt/diffuse-shared/<your-experiment-dir>
md_workflows.mdmx \
  --param-pdb-id 6B8X \
  --ix 1 \
  --ntomp 26 \
  --resolv-ntmpi 8 \
  --resolv-ntomp 1
```

To build the ACTL overlay locally:

```bash
docker buildx build --platform linux/amd64 \
  -f Dockerfile.astera \
  --build-arg MD_WORKFLOWS_BASE_IMAGE=docker.io/diffuseproject/md:0.0.2@sha256:0ec5455d36f3d097fa67c73c3b7b86c0bd039ec19d5ca416ae5350d9476703e5 \
  -t harbor.astera.sh/library/md-workflows:0.0.2-actl-2026-06-08 \
  .
```

## 1) Build the Docker image

From the project root (where the `Dockerfile` is):

```bash
docker pull diffuseproject/md:0.0.2
```

## 2) Start a container

Run the container interactively, mounting the current project directory so inputs/outputs are available on your host:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  --gpus all \
  --name md_container \
  -e HOME=/workspace \
  -v "$(pwd):/workspace" \
  -w /workspace \
  diffuseproject/md:0.0.2 \
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
