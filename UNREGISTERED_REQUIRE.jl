PKGS_USED_CLONABLE = ["https://github.com/seung-lab/ImageRegistration.git", 
        		      "https://github.com/madeleineudell/ParallelSparseMatMul.jl.git",
                      "https://github.com/seung-lab/CloudVolume.jl.git"
                      ]

for pkg in PKGS_USED_CLONABLE
	Pkg.clone(pkg)
end