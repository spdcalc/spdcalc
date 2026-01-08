#!/usr/bin/env bash

set -euxo pipefail
echo "Setting up development environment..."

# Persistent Zsh history setup https://gavinest.com/posts/devcontainers-persist-zsh-history/
sudo mkdir -p /commandhistory
sudo touch /commandhistory/.zsh_history
sudo chown -R ${USER} /commandhistory

echo "autoload -Uz add-zsh-hook; append_history() { fc -W }; add-zsh-hook precmd append_history; export HISTFILE=/commandhistory/.zsh_history" >> /home/${USER}/.zshrc
