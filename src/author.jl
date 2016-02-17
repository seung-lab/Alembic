function author()
	return Dict{Any, Any}(
					"by"	  => ENV["USER"],
					"machine" => gethostname(),
					"timestamp" => string(now())
		)
end

