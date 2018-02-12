PKGS_USED_CLONABLE = ["https://github.com/seung-lab/ImageRegistration.git", 
                      "https://github.com/seung-lab/CloudVolume.jl.git",
                      "https://github.com/jingpengw/AWSCore.jl.git",
                      "https://github.com/jingpengw/AWSSDK.jl.git",
                      "https://github.com/jingpengw/AWSSQS.jl.git"
                      ]

for pkg in PKGS_USED_CLONABLE
	Pkg.clone(pkg)
end