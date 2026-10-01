function main_all_benchmarks()
% MAIN_ALL_BENCHMARKS Runs all 7 benchmark problems from Ko, Lee, Bathe (2017)
% "The MITC4+ shell element in geometric nonlinear analysis"

fprintf('=================================================================\n');
fprintf('  MITC4+ SHELL ELEMENT IN GEOMETRIC NONLINEAR ANALYSIS (2017)\n');
fprintf('  MATLAB Benchmark Suite - All 7 Numerical Examples\n');
fprintf('=================================================================\n\n');

try
    run_cantilever();
    run_slit_annular();
    run_hemisphere();
    run_scordelis();
    run_cylinder();
    run_twisted_cantilever();
    run_hook();

    fprintf('=================================================================\n');
    fprintf('  ALL 7 BENCHMARK PROBLEMS COMPLETED SUCCESSFULLY!\n');
    fprintf('=================================================================\n');
catch ME
    fprintf('Error running benchmarks: %s\n', ME.message);
    rethrow(ME);
end

end
