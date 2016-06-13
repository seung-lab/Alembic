function sync_images(src_image_ref, dst_image_ref)
	src_image_local = fetch(src_image_ref);
	dst_image_local = fetch(dst_image_ref);
	global SHARED_SRC_IMAGE = SharedArray(eltype(src_image_local), size(src_image_local), pids=local_procs());
	global SHARED_DST_IMAGE = SharedArray(eltype(dst_image_local), size(dst_image_local), pids=local_procs());
	SHARED_SRC_IMAGE[:, :] = src_image_local[:, :];
	SHARED_DST_IMAGE[:, :] = dst_image_local[:, :];

	for pid in local_procs()
	remotecall(pid, sync_images_subroutine, SHARED_SRC_IMAGE, SHARED_DST_IMAGE);
      end
#=	tofetch = Array{RemoteRef}(0);
	for pid in local_procs()
		push!(tofetch, remotecall(pid, sync_images_subroutine, SHARED_SRC_IMAGE, SHARED_DST_IMAGE));
	end
      	for ref in tofetch
		fetch(ref);
	end
=#
end

function sync_images_subroutine(local_src_image, local_dst_image)
	global SHARED_SRC_IMAGE = local_src_image;
	global SHARED_DST_IMAGE = local_dst_image;
end

function my_host_addr()
	return Base.Worker(myid()).bind_addr;
end

function get_local_host_addr(id)
	return remotecall_fetch(id, my_host_addr);
end

function local_procs()
	localhost = Base.Worker(myid()).bind_addr;
	remotehosts = map(get_local_host_addr, procs());
	local_procs_indices = find(p -> p == localhost, remotehosts);
	return procs()[local_procs_indices];
end

function Base.size(r::RemoteRef, args...)
      if r.where == myid()
	return size(fetch(r), args...)
		    end
	return remotecall_fetch(r.where, size, r, args...)
end

function Base.eltype(r::RemoteRef, args...)
      if r.where == myid()
	return eltype(fetch(r), args...)
		    end
	return remotecall_fetch(r.where, eltype, r, args...)
      end

function get_proc_range(idx, arr::Array)
        worker_procs = setdiff(procs(), myid());
	nchunks = length(worker_procs);
if nchunks == 0 return 1:length(arr); end
if idx == myid() return 1:0; end
	splits = [round(Int64, s) for s in linspace(0, length(arr), nchunks + 1)];
	return splits[findfirst(worker_procs, idx)]+1:splits[findfirst(worker_procs, idx) + 1]
	end
