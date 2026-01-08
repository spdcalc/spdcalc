#!/usr/bin/env bash
# Initialize command - runs on HOST before container creation
# Ensures required files exist to prevent mount failures

set -e

# Warn if SSH agent is not available
# if [ -z "$SSH_AUTH_SOCK" ]; then
#     echo "Warning: SSH_AUTH_SOCK not set. SSH agent forwarding will not work."
#     echo "Start ssh-agent on your host:"
#     echo "  eval \"\$(ssh-agent -s)\""
#     echo "  ssh-add ~/.ssh/id_ed25519"
# fi
