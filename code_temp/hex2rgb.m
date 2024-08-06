% (c) Jakob Nikolas Kather 2015
% contact: http://www.kather.me
% https://gist.github.com/jnkather/6de1287c446713266e63

function rgb = hex2rgb(hexString)
	if size(hexString,2) ~= 6
		error('invalid input: not 6 characters');
	else
		r = double(hex2dec(hexString(1:2)))/255;
		g = double(hex2dec(hexString(3:4)))/255;
		b = double(hex2dec(hexString(5:6)))/255;
		rgb = [r, g, b];
	end
end