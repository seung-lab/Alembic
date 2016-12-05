using Base.Test

image_without_mask_fn = "test_images/img_without_mask_test.h5"
image_with_mask_fn = "test_images/img_with_mask_test.h5"

@test image_has_mask(image_without_mask_fn) == false
@test image_has_mask(image_with_mask_fn) == true

image_without_mask = get_image(image_without_mask_fn)
image_with_mask = get_image(image_with_mask_fn)

@test image_has_mask(image_without_mask) == false
@test image_has_mask(image_with_mask) == true

mask_without = get_image_mask(image_without_mask_fn)
mask = get_image_mask(image_with_mask_fn)

@test mask_without .== trues(size(image_without_mask))