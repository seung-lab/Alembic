using Alembic
using DataFrames
image_dim = Int64(30000/4)

pt = (0,0)
ranges = (199, (1:227, 2150:2550), [27, 201], 200, (1:527, 1850:2850), (-473:527, 1850:2850), [27, 501], [501, 501], [0, 0])

src_img = SharedArray(rand(UInt8, image_dim, image_dim))
dst_img = SharedArray(rand(UInt8, image_dim, image_dim))
dst_image = dst_img
src_image = src_img
#(199, (1:227, 2150:2550), [27, 201], 200, (1:527, 1850:2850), (-473:527, 1850:2850), [27, 501], [501, 501], [0, 0])
#(199, (1:227, 2350:2750), [27, 201], 200, (1:527, 2050:3050), (-473:527, 2050:3050), [27, 501], [501, 501], [0, 0])
#(199, (1:227, 2550:2950), [27, 201], 200, (1:527, 2250:3250), (-473:527, 2250:3250), [27, 501], [501, 501], [0, 0])
#(199, (1:227, 2750:3150), [27, 201], 200, (1:527, 2450:3450), (-473:527, 2450:3450), [27, 501], [501, 501], [0, 0])
#(199, (1:227, 2950:3350), [27, 201], 200, (1:527, 2650:3j650), (-473:527, 2650:3650), [27, 501], [501, 501], [0, 0])
#(199, (1:227, 3150:3550), [27, 201], 200, (1:527, 2850:3850), (-473:527, 2850:3850), [27, 501], [501, 501], [0, 0])
#(199, (1:227, 3350:3750), [27, 201], 200, (1:527, 3050:4050), (-473:527, 3050:4050), [27, 501], [501, 501], [0, 0])

A = SharedArray(rand(UInt8, 1001, 1001))
B = SharedArray(rand(UInt8, 227, 401))

N=10
src_index, src_range, src_pt_loc, dst_index, dst_range, dst_range_full, dst_pt_loc, dst_pt_loc_full, rel_offset = ranges;

dst_quart_range_i = linspace(dst_range[1][1], dst_range[1][end], 5)
dst_quart_range_j = linspace(dst_range[2][1], dst_range[2][end], 5)
cent_sum = 0;
cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]):round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[2])]);
cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]):round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[4])]);
cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[2]), round(Int64, dst_quart_range_j[2]):round(Int64, dst_quart_range_j[4])]);
cent_sum += sum(dst_image[round(Int64, dst_quart_range_i[4]), round(Int64, dst_quart_range_j[2]):round(Int64, dst_quart_range_j[4])]);

correspondence_properties = DataFrame();
correspondence_properties[:ranges_src_pt_loc] = [src_pt_loc];
correspondence_properties[:ranges_src_range] = src_range;
correspondence_properties[:ranges_dst_pt_loc] = [dst_pt_loc];
correspondence_properties[:ranges_dst_range] = dst_range;
correspondence_properties[:ranges_dst_range_full] = dst_range_full;
#correspondence_properties[:ranges_scale] = scale;

src_patch, dst_patch = prepare_patches(src_image, dst_image, src_range, dst_range, dst_range_full)


@time begin
for i in 0:100
   # get_match(pt, ranges, src_img, dst_img)
   #convolve_Float64_planned(A, B)
   normxcorr2_preallocated(A, B)
end
end
println("Done")
