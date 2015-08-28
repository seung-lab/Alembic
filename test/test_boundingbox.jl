using Base.Test

# test_bb_operations()
A = BoundingBox(10,10,30,20)
B = BoundingBox(5,20,20,10)
C = A+B
D = BoundingBox(5,10,35,20)
@test C == D

@test_throws InexactError BoundingBox(0,0,-1,-1)

# test_find_mesh_bb()
  dst = [4.0 2.0;
          8.0 2.0;
          6.0 10.0]
  bounds = find_bounds(dst')
  @test bounds == (3, 1, 9, 11)

  dst = [0.0 0.0;
          0.0 0.0;
          0.0 0.0]
  bounds = find_bounds(dst')
  @test bounds == (-1, -1, 1, 1)

# test_tform_bb()
sz = (100, 100)
tform = [1 0 0;
        0 1 0;
        100 100 1];
bb = BoundingBox(100, 100, 100, 100)
tbb = tform_bb(sz2bb(sz), tform)
@test bb2pts(bb) == bb2pts(tbb)

sz = (100, 100)
tform = [cos(pi/6) -sin(pi/6) 0;
        sin(pi/6) cos(pi/6) 0;
        0 0 1];
bb = BoundingBox(0, -50, 136.60254037844388, 136.60254037844388)
tbb = tform_bb(sz2bb(sz), tform)
@test_approx_eq bb2pts(bb) bb2pts(tbb)

# test_snap_bb()
bb = snap_bb(sz2bb((100.5, 200.75)))
tbb = sz2bb((101, 201))
@test bb2pts(bb) == bb2pts(tbb)

bb = snap_bb(sz2bb((100.9, 200.1)))
tbb = sz2bb((101, 201))
@test bb2pts(bb) == bb2pts(tbb)

bb = snap_bb(BoundingBox(100.9, 200.1, 100.5, 100.2))
tbb = sz2bb((101, 201))
tbb = BoundingBox(100, 200, 102, 101)
@test bb2pts(bb) == bb2pts(tbb)  
end