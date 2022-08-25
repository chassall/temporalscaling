%% helper function for gaussian processes

function sqexp = squared_exponential(distances, length_scale)

sqexp = exp((distances/length_scale).^2/-2);

end