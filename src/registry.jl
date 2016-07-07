function get_registries()
	return REGISTRY_PREMONTAGED, REGISTRY_MONTAGED, REGISTRY_PREALIGNED, REGISTRY_ALIGNED;
end

function sync_registries()
	ref = RemoteRef();
	put!(ref, get_registries())

	procs_to_call = setdiff(procs(), [myid()])
	pmap(remotecall_fetch, procs_to_call, repeated(sync_registries), repeated(ref))

	close(ref);
end

function sync_registries(ref)
    	REGISTRY_PREMONTAGED, REGISTRY_MONTAGED, REGISTRY_PREALIGNED, REGISTRY_ALIGNED = fetch(ref);
end


