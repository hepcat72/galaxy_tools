#!/bin/bash

#    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
#    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#  Modified by Lance Parsons (lparsons@princeton.edu)
#  2011-03-15  Adapted to allow galaxy to determine filetype
#  2015-10-21  Updated to make compatible with OSX (BSD sed)
#  2015-11-13  Removed LIBRARY_NAME, no longer needed
#  2016-04-28  Output summary as simple tabular output

# This is a shell script wrapper for 'fastx_barcode_splitter.pl'
#
# 1. Output files are saved at the dataset's files_path directory.

if [ "$1x" = "x" ]; then
  echo "Usage: $0 [OUTPUT_PATH] [OUTFILETYPE] [OPTIONS]" >&2
  exit 1
fi

OUTPUT_PATH="$1"
FILETYPE="$2"
shift 2
# The rest of the parameters are passed to the split program

if [ "${OUTPUT_PATH}x" = "x" ]; then
  echo "Usage: $0 [OUTPUT_PATH] [OUTFILETYPE] [OPTIONS]" >&2
  exit 1
fi

mkdir -p "$OUTPUT_PATH"
if [ ! -d "$OUTPUT_PATH" ]; then
  echo "Error: failed to create output path '$OUTPUT_PATH'" >&2
  exit 1
fi

PREFIX="$OUTPUT_PATH/"
SUFFIX=".$FILETYPE"
DIRECTORY="$(cd "$(dirname "$0")" && pwd)"

RESULTS=$("$DIRECTORY/barcode_splitter.py" --format "$FILETYPE" --prefix "$PREFIX" --suffix "$SUFFIX" "$@")
if [ $? != 0 ]; then
  echo "error"
fi

echo "$RESULTS"