# MPLAPACK Development Environment with Local LLM

This directory provides Docker configurations for building and running MPLAPACK,
optionally integrated with a fully offline Qwen3.5 LLM coding assistant.

## Files

```
docker/devel/
├── Dockerfile-ubuntu24.04-devel   # MPLAPACK dev image (CUDA, Intel oneAPI, MinGW, aider, Emacs)
├── docker-compose-llm.yml         # Qwen3.5 inference server (llama.cpp + CUDA)
├── .env.llm.example               # Configuration template for the inference server
├── download.py                    # Model downloader script
└── README.md                      # This file
```

---

## 1. MPLAPACK Development Container

### Build

```bash
# Public users (read-only, HTTPS)
docker build -f Dockerfile-ubuntu24.04-devel -t mplapack:devel2404 .

# Maintainers (read-write, SSH)
DOCKER_BUILDKIT=1 docker build \
  -f Dockerfile-ubuntu24.04-devel \
  --build-arg DOCKER_UID="$(id -u)" \
  --build-arg DOCKER_GID="$(id -g)" \
  --build-arg USE_SSH=true \
  --build-arg SSH_KEY_NAME=id_ed25519 \
  --secret id=ssh_key,src="$HOME/.ssh/id_ed25519" \
  -t mplapack:devel2404 .

# Fork users (custom repo, HTTPS)
DOCKER_BUILDKIT=1 docker build \
  -f Dockerfile-ubuntu24.04-devel \
  --build-arg DOCKER_UID="$(id -u)" \
  --build-arg DOCKER_GID="$(id -g)" \
  --build-arg MPLAPACK_REPO=https://github.com/yourname/mplapack.git \
  -t mplapack:devel2404 .
```

### Run

```bash
# Basic shell
docker run --privileged --gpus all -it --rm mplapack:devel2404 /bin/bash

# With SSH agent for git push/pull (maintainers)
eval $(ssh-agent)
ssh-add ~/.ssh/id_ed25519
docker run --privileged --gpus all -it --rm \
  -v $SSH_AUTH_SOCK:/ssh-agent \
  -e SSH_AUTH_SOCK=/ssh-agent \
  mplapack:devel2404 /bin/bash

# With LLM assistant (start inference server first; see section 2)
docker run --privileged --gpus all -it --rm \
  --network mplapack-llm_llm-net \
  -e OPENAI_API_BASE=http://172.31.0.10:8080/v1 \
  -e NO_COLOR=1 \
  mplapack:devel2404 /bin/bash
```

### What is included

| Component | Purpose |
|---|---|
| GCC / GFortran / Clang | C++/Fortran build toolchain |
| Intel oneAPI (icx, ifx, MKL, MPI) | Intel compiler and math libraries |
| MinGW-w64 + Wine | Windows cross-compilation and test execution |
| CUDA 13.1.1 devel | GPU-accelerated MPLAPACK builds |
| aider 0.85.2 | LLM coding agent CLI |
| Emacs (nox) + aider.el | Editor with LLM integration, no X11 required |
| universal-ctags, ripgrep | Repository indexing for aider repo-map |
| Japanese locale (ja_JP.UTF-8) | Japanese input/output in terminal and Emacs |

---

## 2. LLM Inference Server

A local OpenAI-compatible API server running Qwen3.5 via llama.cpp.
After the initial model download, no internet connection is required.

### Prerequisites

- NVIDIA GPU with >= 40GB VRAM (A100 80GB recommended)
- [nvidia-container-toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) installed on the host
- Model file downloaded (see below)

### Step 1: Download the Model (one-time, requires internet)

```bash
# Create a temporary venv for the downloader
python3 -m venv /tmp/hf-venv
/tmp/hf-venv/bin/pip install -q huggingface_hub

# Option A (recommended): 35B-A3B Q8_0 (~37GB) - near-lossless quality, fast
/tmp/hf-venv/bin/python3 download.py --model A --dest /work1/llm_models/Qwen/

# Option B: 35B-A3B Q4_K_M (~20GB) - good quality, fastest
/tmp/hf-venv/bin/python3 download.py --model B --dest /work1/llm_models/Qwen/

# Option C: 122B-A10B UD-Q3_K_XL (~53GB) - best quality, slower  [split files]
/tmp/hf-venv/bin/python3 download.py --model C --dest /work1/llm_models/Qwen/

# Option D: 122B-A10B UD-Q4_K_XL (~69GB) - best quality, larger  [split files]
/tmp/hf-venv/bin/python3 download.py --model D --dest /work1/llm_models/Qwen/

# Discard venv (model file is already saved)
rm -rf /tmp/hf-venv
```

### Model Selection Guide (A100 80GB)

| Option | Model | Quant | File size | VRAM | Quality | Speed |
|---|---|---|---|---|---|---|
| A | 35B-A3B | Q8_0 | ~37GB | ~39GB | near-lossless | fast |
| B | 35B-A3B | Q4_K_M | ~20GB | ~22GB | good | fastest |
| C | 122B-A10B | UD-Q3_K_XL | ~53GB | ~56GB | best | slower |
| D | 122B-A10B | UD-Q4_K_XL | ~69GB | ~72GB | best | slower |

Options C and D use the Unsloth Dynamic 2.0 format and are stored as split GGUF files
in a subdirectory. Set `MODEL_FILE` to the first shard filename
(e.g. `UD-Q3_K_XL/Qwen3.5-122B-A10B-UD-Q3_K_XL-00001-of-XXXX.gguf`).
The exact path is printed by `download.py` at the end of the download as `Load with llama-server:`.

### Step 2: Configure and Start

```bash
cp .env.llm.example .env.llm
# Edit .env.llm: set MODEL_DIR to the directory containing the GGUF file,
# and set MODEL_FILE to match the downloaded filename.
nano .env.llm

docker compose -f docker-compose-llm.yml --env-file .env.llm build inference
docker compose -f docker-compose-llm.yml --env-file .env.llm up -d

# Wait for model load (60-120 seconds)
docker compose -f docker-compose-llm.yml logs -f | grep -m1 "server is listening"

# Verify from host
curl http://localhost:8080/health
# Expected: {"status":"ok"}
```

### Stop

```bash
docker compose -f docker-compose-llm.yml down
```

---

## 3. Using aider in the Dev Container

### CLI

```bash
# Inside the dev container connected to llm-net:
cd ~/mplapack
aider \
  --model openai/Qwen3.5-35B-A3B \
  --openai-api-base http://172.31.0.10:8080/v1 \
  --openai-api-key dummy \
  --no-check-update \
  --no-show-model-warnings \
  --no-pretty \
  --no-fancy-input \
  --map-tokens 1024 \
  --max-chat-history-tokens 4096 \
  --edit-format diff \
  --auto-commits
```

Useful commands inside the aider session:

| Command | Action |
|---|---|
| `/add src/Rlasq2.cpp` | Add file to context |
| `/ask <question>` | Ask without editing files |
| `/diff` | Show diff of last change |
| `/undo` | Revert last commit |
| `/tokens` | Show current token usage |
| `/clear` | Reset chat history |
| `/quit` | Exit |

### Emacs + aider.el

```bash
# Launch Emacs with a source file
docker run --privileged --gpus all -it --rm \
  --network mplapack-llm_llm-net \
  -e NO_COLOR=1 \
  -v /path/to/mplapack:/home/docker/mplapack:rw \
  mplapack:devel2404 emacs ~/mplapack/src/Rlasq2.cpp
```

Key bindings (pre-configured in `~/.emacs.d/init.el`):

| Key | Action |
|---|---|
| `C-c a a` | Start aider session |
| `C-c a f` | Add current file to aider context |
| `C-c a q` | Ask a question (no edit) |
| `C-c a c` | Request a code change |
| `C-c a d` | Show diff of last change |
| `C-c a u` | Undo last aider change |

---

## 4. Offline Operation

After the initial setup (model download + image build), everything runs
without internet access.

```bash
# Save images on an internet-connected machine
docker save mplapack:devel2404 | gzip > mplapack_devel2404.tar.gz
docker save ghcr.io/ggerganov/llama.cpp:server-cuda | gzip > llamacpp_server_cuda.tar.gz

# Transfer to air-gapped machine (rsync / USB / etc.)
rsync -avP mplapack_devel2404.tar.gz user@airgapped:~/
rsync -avP llamacpp_server_cuda.tar.gz user@airgapped:~/
rsync -avP /path/to/models/Qwen3.5-35B-A3B-Q8_0.gguf user@airgapped:/path/to/models/

# Load on air-gapped machine
docker load < mplapack_devel2404.tar.gz
docker load < llamacpp_server_cuda.tar.gz
```

---

## 5. Troubleshooting

**GPU not visible in container:**
```bash
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi
```

**Model architecture error (`unknown model architecture: 'qwen35moe'`):**
The pre-built llama.cpp image may be too old. Build from source using
`Dockerfile.inference` from the `llm-local-agent/` directory and replace
the `image:` line in `docker-compose-llm.yml` with a `build:` section.

**Context window exceeded in aider:**
Increase `LLAMA_CTX_SIZE` in `.env.llm` (A100 80GB supports up to ~500K
for 35B-Q8) and restart the inference server.

**Japanese not displaying correctly:**
Ensure your terminal uses UTF-8 encoding. The container sets
`LANG=ja_JP.UTF-8` and Emacs is configured for UTF-8 by default.
