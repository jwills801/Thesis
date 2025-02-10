{ pkgs, lib, config, inputs, ... }:

{
  # https://devenv.sh/basics/
  env.GREET = "devenv";

  # https://devenv.sh/packages/
  packages = [ pkgs.git  
               pkgs.pre-commit
               pkgs.libGLU           # necessary for wecopttool.geom
               pkgs.libGL            # necessary for wecopttool.geom
               pkgs.xorg.libXrender  # necessary for wecopttool.geom
               pkgs.xorg.libXcursor  # necessary for wecopttool.geom
               pkgs.xorg.libXfixes   # necessary for wecopttool.geom
               pkgs.xorg.libXft      # necessary for wecopttool.geom
               pkgs.fontconfig       # necessary for wecopttool.geom
               pkgs.xorg.libXinerama # necessary for wecopttool.geom
               pkgs.xorg.libX11      # necessary for wecopttool.geom
               ];          

# configure pre-commit to do the things we want (format, lint, and clear)
git-hooks = {
  hooks = {
    ruff.enable = true; # linter and formater for python
  };

  # Configure nbstripout to run everytime we commit a jupyter notebook file
  hooks.nbstripout = {
    enable = true;

    # The name of the hook (appears on the report table):
    name = "nbstripout";

    # The command to execute (mandatory):
    entry = "nbstripout";

    # The pattern of files to run on
    files = "\\.ipynb$";

    # List of file types to run on
    types = [ "text" ];

    language = "python";
};

  # https://devenv.sh/languages/
    # Enable Python
    languages.python.enable = true;

    # Create a virtual environment and install wecopttool
    languages.python.venv = {
      enable = true;
      requirements = ''
      wecopttool[geometry]
      nbstripout # strips my jupyter notebooks of outputs before committing
      jupyter
      '';
    };

  # https://devenv.sh/processes/
  # processes.cargo-watch.exec = "cargo-watch";

  # https://devenv.sh/services/
  # services.postgres.enable = true;

  # https://devenv.sh/scripts/
  scripts.hello.exec = ''
    echo hello from $GREET
  '';

  enterShell = ''
    hello
    git --version
  '';

  # https://devenv.sh/tasks/
  # tasks = {
  #   "myproj:setup".exec = "mytool build";
  #   "devenv:enterShell".after = [ "myproj:setup" ];
  # };

  # https://devenv.sh/tests/
  enterTest = ''
    echo "Running tests"
    git --version | grep --color=auto "${pkgs.git.version}"
  '';

  # https://devenv.sh/pre-commit-hooks/
  # pre-commit.hooks.shellcheck.enable = true;

  # See full reference at https://devenv.sh/reference/options/
}
