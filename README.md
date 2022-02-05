# chartools-dev
Suite of tools to work with ChAR-seq data. 

## Installation

### Install the python chartools package
Clone this repository then `pip install <PATH_TO_REPOSITORY_CLONE>`


### Install the Julia package
Install Julia 1.6 or 1.7 : https://julialang.org/downloads/, then add Julia to search path `export PATH="$PATH:/path/to/<Julia directory>/bin"`

Then we need to install a few Julia packages. First launch the Julia interactive shell

`julia`

Then within julia, press `]` to enter pkg mode, and install the necessary packages by entering
- `add FASTX` (https://github.com/BioJulia/FASTX.jl)
- `add CodecZlib` (https://github.com/JuliaIO/CodecZlib.jl)
- `add TranscodingStreams` (https://github.com/JuliaIO/TranscodingStreams.jl) 


The first package is part of BioJulia to read fastq files, and the second two packages are to read/write gzipped files. You can now quit the interactive sheel with `CTRL-D` or `quit()`

Finally create a file in your home directory `~/.juliarc.jl`, with the following content, to allow Julia to find this Jchartools package

```
push!(LOAD_PATH,"/path/to/this_repo_clone/Jchartools")
```

