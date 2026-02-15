#!/usr/bin/env bash
# Make a small preview GIF from sample_video.webm so GitHub can display it inline.
# Requires ffmpeg. Usage: ./scripts/make_preview_gif.sh [input.webm]
# Output: sample_video_preview.gif (first 6 s, 480px wide, 8 fps, palette-optimized)

set -e
INPUT="${1:-sample_video.webm}"
OUTPUT="${INPUT%.*}_preview.gif"

if ! command -v ffmpeg &>/dev/null; then
  echo "Need ffmpeg: brew install ffmpeg"
  exit 1
fi
if [[ ! -f "$INPUT" ]]; then
  echo "Input not found: $INPUT"
  exit 1
fi

echo "Creating palette..."
ffmpeg -y -i "$INPUT" -vf "fps=8,scale=480:-1:flags=lanczos,palettegen" palette.png

echo "Encoding GIF (first 6 s)..."
ffmpeg -y -i "$INPUT" -i palette.png -lavfi "fps=8,scale=480:-1:flags=lanczos[x];[x][1:v]paletteuse" -t 6 "$OUTPUT"

rm -f palette.png
echo "Done: $OUTPUT"
