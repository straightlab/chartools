# chartools-dev
Suite of tools to work with ChAR-seq data. 

## Installation

### Install the python chartools package
Clone this repository then `pip install <PATH_TO_REPOSITORY_CLONE>`


### Install the Julia package
This is a mess right now, as this package was written for an older version of Julia, v0.6.

To install Julia 0.6, the easiest way is by downloading the precompiled binary, as described here
`https://julialang.org/downloads/oldreleases/`
and here
`https://julialang.org/downloads/platform/#linux_and_freebsd`

Basically, for Linux this boils down to 
- downloading the binary `wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz`
- unpack `tar zxvf julia-0.6.4-linux-x86_64.tar.gz`
- add Julia to search path `export PATH="$PATH:/path/to/<Julia directory>/bin"`

Then we need to install a few Julia packages. First launch the Julia interactive shell

`julia`

Then within julia run:
- `Pkg.add("BioSequences")`
- `Pkg.add("CodecZlib")` 
- `Pkg.add("TranscodingStreams")`

The first package is part of BioJulia to read fastq files, and the second two packages are to read/write gzipped files. You can now quit the interactive sheel with `CTRL-D` or `quit()`

Finally create a file in your home directory `~/.juliarc.jl`, with the following content

```
push!(LOAD_PATH,"/path/to/this_repo_clone/Jchartools")
```

