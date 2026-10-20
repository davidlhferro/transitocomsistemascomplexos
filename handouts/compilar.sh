#!/usr/bin/env bash
# Compila todos os handouts e coloca os PDFs em handouts/pdf/
set -e
cd "$(dirname "$0")"
mkdir -p pdf
export TEXINPUTS="$PWD/estilo//:$PWD/../figuras//:"
for f in teoria/*.tex projeto/*.tex; do
  [ -e "$f" ] || continue
  latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build "$f" > /dev/null
  cp "build/$(basename "${f%.tex}").pdf" pdf/
  echo "ok: $f"
done
