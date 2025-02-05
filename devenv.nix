{ pkgs, lib, config, inputs, ... }:

{
  # https://devenv.sh/basics/
  env.GREET = "devenv";

  # https://devenv.sh/packages/
  packages = [ pkgs.git  
               pkgs.pre-commit
               ];          

  # https://devenv.sh/languages/
    # Enable Python
    languages.python.enable = true;

    # Create a virtual environment and install wecopttool
    languages.python.venv = {
      enable = true;
      requirements = ''
        wecopttool
        nbstripout
      '';
    };
    # nbstripout is for removing the output of jupyter notebooks before they are commited
      # this saves space in git and makes it so you have to regenerate the output - saving potential errors
      # The nbstripout is made a hook in the file named .pre-commit-config.yaml

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
