function [psnr,ssim,sam] = quality_assessment(rec,orig)
rec  = double(im2uint8(rec));
orig = double(im2uint8(orig));
psnr = psnr_ours(rec,orig);
ssim = ssim_ours(rec,orig);
sam  = sam_ours(rec,orig);
end
