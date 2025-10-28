{ pkgs ? import (builtins.fetchTarball {
    name = "nixpkgs_24.11_2025-05-05";
    url = "https://github.com/NixOS/nixpkgs/archive/5b35d248e9206c1f3baf8de6a7683fee126364aa.tar.gz";
    sha256 = "11sih5m3cqmc9dfsjqj51pn3f26bjbhx1d0pgm7sggyh68wllfrm";
  }) {} }:

pkgs.mkShell {
  nativeBuildInputs = with pkgs.buildPackages; [
    gcc
    zlib
    python3
  ];
}
