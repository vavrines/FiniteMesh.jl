using PyCall, Conda

# try importing meshio
# catch the installation
# 1) Julia built-in miniconda
# 2) global pip installer

#Conda.add_channel("conda-forge")
#Conda.add("meshio")
#Conda.add("scipy")
#Conda.pip_interop(true)
#Conda.pip("install", "meshio")

@info "installing python mates"
cmd = `pip3 install meshio scipy --user`
run(cmd)

#meshio = pyimport("meshio")
#scipy = pyimport("scipy")
