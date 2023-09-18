#!/bin/bash
# put current commit into c
#$cur_commit=git rev-parse HEAD
# $status=git status -s
cat > AeF-hyperfine-structure.inl <<EOF
// AeF-hyperfine-structure inline props
// this is an auto-generated file, do not modify
constexpr char aef_git_commit[]="$(git rev-parse HEAD)";
constexpr char aef_git_status[]=R"/RAWCHARS($(git status -s))/RAWCHARS";
EOF
