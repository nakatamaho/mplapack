#!/usr/bin/env python3
# download.py
# Download Qwen3.5 GGUF model files from Hugging Face.
#
# Usage:
#   python3 -m venv /tmp/hf-venv
#   /tmp/hf-venv/bin/pip install -q huggingface_hub
#   /tmp/hf-venv/bin/python3 download.py --model A --dest /path/to/models
#   rm -rf /tmp/hf-venv
#
# Notes:
#   35B models   : single GGUF file, downloaded with hf_hub_download
#   122B models  : split GGUF files in a subdirectory (Unsloth Dynamic 2.0),
#                  downloaded with snapshot_download + allow_patterns

import argparse
import sys
from pathlib import Path

MODELS = {
    "A": {
        "desc": "35B-A3B Q8_0  (~37GB) - near-lossless quality, fast",
        "repo_id": "unsloth/Qwen3.5-35B-A3B-GGUF",
        "filename": "Qwen3.5-35B-A3B-Q8_0.gguf",
        "split": False,
    },
    "B": {
        "desc": "35B-A3B Q4_K_M (~20GB) - good quality, fastest",
        "repo_id": "unsloth/Qwen3.5-35B-A3B-GGUF",
        "filename": "Qwen3.5-35B-A3B-Q4_K_M.gguf",
        "split": False,
    },
    "C": {
        "desc": "122B-A10B UD-Q3_K_XL (~53GB) - best quality, slower  [split files]",
        "repo_id": "unsloth/Qwen3.5-122B-A10B-GGUF",
        "pattern": "UD-Q3_K_XL/*",
        "split": True,
    },
    "D": {
        "desc": "122B-A10B UD-Q4_K_XL (~69GB) - best quality, larger  [split files]",
        "repo_id": "unsloth/Qwen3.5-122B-A10B-GGUF",
        "pattern": "UD-Q4_K_XL/*",
        "split": True,
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download Qwen3.5 GGUF model from Hugging Face.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(f"  {k}: {v['desc']}" for k, v in MODELS.items()),
    )
    parser.add_argument(
        "--model",
        choices=list(MODELS.keys()),
        default="A",
        help="Model variant to download (default: A)",
    )
    parser.add_argument(
        "--dest",
        required=True,
        type=Path,
        help="Local directory to save the model file(s).",
    )
    args = parser.parse_args()

    try:
        from huggingface_hub import hf_hub_download, snapshot_download
    except ImportError:
        print("ERROR: huggingface_hub is not installed.", file=sys.stderr)
        print("  pip install huggingface_hub", file=sys.stderr)
        sys.exit(1)

    model = MODELS[args.model]
    dest = args.dest
    dest.mkdir(parents=True, exist_ok=True)

    print(f"Model   : {model['desc']}")
    print(f"Repo    : {model['repo_id']}")
    print(f"Dest    : {dest}")
    print("Downloading... (resumable)\n")

    if model["split"]:
        # 122B: split GGUF files in subdirectory - use snapshot_download
        print(f"Pattern : {model['pattern']}")
        path = snapshot_download(
            repo_id=model["repo_id"],
            allow_patterns=[model["pattern"]],
            local_dir=str(dest),
        )
        total_gb = sum(
            f.stat().st_size for f in Path(path).rglob("*.gguf")
        ) / 1024**3
        subdir = model["pattern"].replace("/*", "")
        print(f"\nDone: {dest}/{subdir}/  (total {total_gb:.1f} GB)")
        print(f"Load with llama-server: --model {dest}/{subdir}/Qwen3.5-122B-A10B-{subdir}-00001-of-*.gguf")
    else:
        # 35B: single GGUF file
        print(f"File    : {model['filename']}")
        path = hf_hub_download(
            repo_id=model["repo_id"],
            filename=model["filename"],
            local_dir=str(dest),
        )
        size_gb = Path(path).stat().st_size / 1024**3
        print(f"\nDone: {path}  ({size_gb:.1f} GB)")
        print(f"Load with llama-server: --model {path}")


if __name__ == "__main__":
    main()
