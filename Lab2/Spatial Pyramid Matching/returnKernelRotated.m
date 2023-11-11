% Function to return rotated Sobel kernel for a given input number of 45deg
% rotations from that returned by default (fpsecial('sobel'))

% Put here to be reused over extracting weak and strong features

function kernel = returnKernelRotated(n_rot)

% Giant, inelegant if-else to specify edge detection kernels
if n_rot == 1
    kernel = fspecial('sobel');
elseif n_rot == 2
    kernel = [0 1 2; -1 0 1; -2 -1 0];
elseif n_rot == 3
    kernel = rot90(fspecial('sobel'), 3);
elseif n_rot == 4
    kernel = flipud([0 1 2; -1 0 1; -2 -1 0]);
elseif n_rot == 5
    kernel = flipud(fspecial('sobel'));
elseif n_rot == 6
    kernel = [0 1 2; -1 0 1; -2 -1 0]';
elseif n_rot == 7
    kernel = fspecial('sobel')';
elseif n_rot == 8
    kernel = flipud([0 1 2; -1 0 1; -2 -1 0]');
end

end

