import numpy as np
import bempp.api

def juffer(dirichl_space, neumann_space, q, x_q, ep_in, ep_ex, kappa):
    from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

    phi_id = sparse.identity(dirichl_space, dirichl_space, dirichl_space)
    dph_id = sparse.identity(neumann_space, neumann_space, neumann_space)
    ep = ep_ex/ep_in

    dF = laplace.double_layer(dirichl_space, dirichl_space, dirichl_space)
    dP = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dirichl_space, kappa)
    L1 = (ep*dP) - dF

    F = laplace.single_layer(neumann_space, dirichl_space, dirichl_space)
    P = modified_helmholtz.single_layer(neumann_space, dirichl_space, dirichl_space, kappa)
    L2 = F - P

    ddF = laplace.hypersingular(dirichl_space, neumann_space, neumann_space)
    ddP = modified_helmholtz.hypersingular(dirichl_space, neumann_space, neumann_space, kappa)
    L3 = ddP - ddF

    dF0 = laplace.adjoint_double_layer(neumann_space, neumann_space, neumann_space)
    dP0 = modified_helmholtz.adjoint_double_layer(neumann_space, neumann_space, neumann_space, kappa)
    L4 = dF0 - ((1.0/ep)*dP0)

    A = bempp.api.BlockedOperator(2, 2)
    A[0, 0] = (0.5*(1.0 + ep)*phi_id) - L1
    A[0, 1] = (-1.0)*L2
    A[1, 0] = L3    # Cambio de signo por definicion de bempp
    A[1, 1] = (0.5*(1.0 + (1.0/ep))*dph_id) - L4

    @bempp.api.real_callable
    def d_green_func(x, n, domain_index, result):
        nrm = np.sqrt((x[0]-x_q[:,0])**2 + (x[1]-x_q[:,1])**2 + (x[2]-x_q[:,2])**2)
        
        const = -1./(4.*np.pi*ep_in)
        result[:] = const*np.sum(q*np.dot(x-x_q, n)/(nrm**3))

    @bempp.api.real_callable
    def green_func(x, n, domain_index, result):
        nrm = np.sqrt((x[0]-x_q[:,0])**2 + (x[1]-x_q[:,1])**2 + (x[2]-x_q[:,2])**2)
        
        result[:] = np.sum(q/nrm)/(4.*np.pi*ep_in)

        
    rhs_1 = bempp.api.GridFunction(dirichl_space, fun=green_func)
    rhs_2 = bempp.api.GridFunction(dirichl_space, fun=d_green_func)

    return A, rhs_1, rhs_2