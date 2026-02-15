#!/usr/bin/env bash
# Convert full video (MP4/WebM) to GIF — whole video, reduced fps to keep size down.
# Requires ffmpeg. Usage: ./scripts/make_preview_gif.sh [input.mp4]
# Output: sample_video.gif (full length, 5 fps, 480px wide, palette-optimized)

set -e
INPUT="${1:-sample_video.mp4}"
OUTPUT="${INPUT%.*}.gif"
FPS=5
WIDTH=480

if ! command -v ffmpeg &>/dev/null; then
  echo "Need ffmpeg: brew install ffmpeg"
  exit 1
fi
if [[ ! -f "$INPUT" ]]; then
  echo "Input not found: $INPUT"
  exit 1
fi

echo "Creating palette (full video, ${FPS} fps, ${WIDTH}px)..."
ffmpeg -y -i "$INPUT" -vf "fps=${FPS},scale=${WIDTH}:-1:flags=lanczos,palettegen" palette.png

echo "Encoding full GIF..."
ffmpeg -y -i "$INPUT" -i palette.png -lavfi "fps=${FPS},scale=${WIDTH}:-1:flags=lanczos[x];[x][1:v]paletteuse" "$OUTPUT"

rm -f palette.png
echo "Done: $OUTPUT"
