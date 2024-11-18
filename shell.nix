{ pkgs ? import <nixpkgs> {} }:
let
  compile-latex = pkgs.writeShellScriptBin "compile-latex" ''
    pdflatex -interaction=nonstopmode reporte.tex
    bibtex reporte
    # pdflatex -interaction=nonstopmode reporte.tex
    pdflatex -interaction=nonstopmode reporte.tex
  '';

  watch-latex = pkgs.writeShellScriptBin "watch-latex" ''
    echo "Watching reporte.tex for changes..."
    ${pkgs.entr}/bin/entr -s 'compile-latex' <<< reporte.tex
  '';
in pkgs.mkShell {
  buildInputs = with pkgs; [
    python3
    python3Packages.pip
    entr
    compile-latex
    watch-latex
  ];

  shellHook = ''
    # Create a virtual environment if it doesn't exist
    if [ ! -d "venv" ]; then
      python3 -m venv venv
    fi

    # Activate the virtual environment
    source venv/bin/activate

    # Upgrade pip
    pip install --upgrade pip

    # Install requirements if requirements.txt exists
    if [ -f requirements.txt ]; then
      pip install -r requirements.txt
    else
      echo "requirements.txt not found. Please create one with your project dependencies."
    fi

    # Inform the user that the environment is ready
    echo "Python virtual environment is set up and requirements are installed."
    echo "To deactivate the virtual environment, type 'deactivate'."
  '';
}