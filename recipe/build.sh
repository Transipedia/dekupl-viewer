#!/bin/bash

mkdir -p "$PREFIX"/usr/local/dekupl-viewer
cp -r src/* "$PREFIX"/usr/local/dekupl-viewer
cp -r toy "$PREFIX"/usr/local

echo "### Install bin ###"
ln -svf "$PREFIX"/usr/local/dekupl-viewer/dekupl-viewer.sh "$PREFIX"/bin/dekupl-viewer
