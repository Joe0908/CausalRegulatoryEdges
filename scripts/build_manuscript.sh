#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_md="$repo_dir/manuscript/When_Does_a_Regulatory_Edge_Become_Causal_initial.md"
output_docx="$repo_dir/manuscript/When_Does_a_Regulatory_Edge_Become_Causal_initial.docx"
output_dir="$repo_dir/manuscript"

if ! command -v pandoc >/dev/null 2>&1; then
  echo "pandoc is required to build the DOCX review copy." >&2
  exit 1
fi

pandoc "$source_md" \
  --from gfm+tex_math_dollars \
  --to docx \
  --resource-path="$repo_dir/manuscript:$repo_dir/reports:$repo_dir" \
  --output "$output_docx"

echo "Wrote $output_docx"

if command -v soffice >/dev/null 2>&1; then
  office_profile="$(mktemp -d)"
  pdf_dir="$(mktemp -d)"
  trap 'rm -rf -- "$office_profile" "$pdf_dir"' EXIT
  soffice -env:UserInstallation="file://$office_profile" --headless \
    --convert-to pdf --outdir "$pdf_dir" "$output_docx" >/dev/null
  mv -f "$pdf_dir/$(basename "${output_docx%.docx}.pdf")" \
    "$output_dir/$(basename "${output_docx%.docx}.pdf")"
  echo "Wrote ${output_docx%.docx}.pdf"
else
  echo "LibreOffice not found; skipped the optional PDF review copy." >&2
fi
