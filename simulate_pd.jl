using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!, u_modified!


function simulate_pd(pd,powergrid,operationpoint,timespan)
    x0 = operationpoint;
    problem = ODEProblem{true}(rhs(powergrid), x0.vec, timespan)
    integrator = init(problem, Rodas4(autodiff=false),saveat=1e-4)

    step!(integrator, pd.tspan_fault[1], true)

    # update integrator with error
    integrator.f = rhs(pd(powergrid))
    u_modified!(integrator,true)

    step!(integrator, pd.tspan_fault[2]-pd.tspan_fault[1], true)

    # update integrator, clear error
    integrator.f = rhs(powergrid)
    u_modified!(integrator,true)

    step!(integrator, timespan[2]-pd.tspan_fault[2], true)


    solve!(integrator)
    result_pd2 = PowerGridSolution(integrator.sol, powergrid)
    return (integrator.sol,result_pd2)
end

function simulate_linefault(disturbed_nodes,powergrid,operationpoint,simulation_time)
	line_fault=LineFault(from=disturbed_nodes[1],to=disturbed_nodes[2]);
	faulty_grid = line_fault(powergrid);
	problem = ODEProblem{true}(rhs(faulty_grid),operationpoint.vec,simulation_time)
	sol_lf = solve(problem,Rodas4(autodiff=false),saveat=1e-4);
	result_lf=PowerGridSolution(sol_lf, faulty_grid);
   	return (sol_lf,result_lf,faulty_grid)
end
