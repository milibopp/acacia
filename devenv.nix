{ pkgs, ... }:

{
  packages = with pkgs; [
    entr
    fd
    git
  ];

  enterShell = ''
  '';

  languages.rust.enable = true;
}
