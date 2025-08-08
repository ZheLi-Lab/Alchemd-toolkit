#!/bin/bash

# Define the patch file and target file
PATCH_FILE="ti.F90.patch"
TARGET_FILE="ti.F90"

# Check if the patch file exists
if [ ! -f "$PATCH_FILE" ]; then
    echo "[Error] Patch file $PATCH_FILE not found. Exiting."
    exit 1
fi

# Check if the target file exists
if [ ! -f "$TARGET_FILE" ]; then
    echo "[Error] Target file $TARGET_FILE not found. Exiting."
    exit 1
fi

# Check if the patch has already been applied by looking for changes in the target file
PATCH_ALREADY_APPLIED=$(patch --dry-run -p0 < "$PATCH_FILE" 2>&1 | grep -i "Reversed (or previously applied)" | wc -l)

if [ "$PATCH_ALREADY_APPLIED" -ne 0 ]; then
    echo "[Info] Patch $PATCH_FILE has already been applied. Skipping."
else
    echo "[Action] Applying patch $PATCH_FILE to $TARGET_FILE..."
    patch -p0 < "$PATCH_FILE"
    if [ $? -ne 0 ]; then
        echo "[Error] Failed to apply patch. Exiting."
        exit 1
    fi
    echo "[Done] Patch applied successfully."
fi

