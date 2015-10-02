for i = 1 : 256^2

	mf(i) = 0;
	c = 0;

	% neigh = zeros(1,4);
	% w = neigh;
	% v = 1e6;

	if i < 256
		if prior.weightNoise(i + 1) < .5			
			mf(i) = mf(i) + prior.av_mess(i + 1);
			% neigh(1) = prior.av_mess(i + 1);
			% w(1) = exp(-(prior.av_mess(i) - prior.av_mess(i + 1) ).^2 ./ (2 * v) );
			c = c + 1;
		end
	end
	if i > 1
		if prior.weightNoise(i - 1) < .5
			mf(i) = mf(i) + prior.av_mess(i - 1);
			% neigh(2) = prior.av_mess(i - 1);
			% w(2) = exp(-(prior.av_mess(i) - prior.av_mess(i - 1) ).^2 ./ (2 * v) );
			c = c + 1;
		end
	end 
	if i > 256
		if prior.weightNoise(i - 256) < .5
			mf(i) = mf(i) + prior.av_mess(i - 256);
			% neigh(3) = prior.av_mess(i - 256);
			% w(3) = exp(-(prior.av_mess(i) - prior.av_mess(i - 256) ).^2 ./ (2 * v) );
			c = c + 1;
		end
	end
	if i < 256 * 255 + 1
		if prior.weightNoise(i + 256) < .5
			mf(i) = mf(i) + prior.av_mess(i + 256);
			% neigh(4) = prior.av_mess(i + 256);
			% w(4) = exp(-(prior.av_mess(i) - prior.av_mess(i + 256) ).^2 ./ (2 * v) );
			c = c + 1;
		end
	end

	if c > 0 
		mf(i) = mf(i) / c;				
		% mf(i) = sum(neigh .* w) ./ sum(w);
	end

end
