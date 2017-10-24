function author()
	return Dict{Any, Any}(
#					"by"	  => ENV["USER"],
					"machine" => gethostname(),
					"timestamp" => string(now())
		)
end

function null_author()
	return Dict{Any, Any}(
					"by"	  => "null",
					"machine" => "null",
					"timestamp" => "null"
		)
end

