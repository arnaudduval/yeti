"""
.. IGA methods
.. by Joaquin Cornejo
.. Disclaimer :: This module will hardly ever be used and 
..               and maybe it does not work as expected
"""

# Python libraries
import numpy as np
from scipy import sparse as sp
import time

# My libraries
from .create_model import thermoMechaModel, eval_basis

def gaussTable(pgaus):
    " Computes Gauss weights and positions in isogeometric space for a known degree "

    if pgaus == 1:
        x = [0.0]
        w = [2.0]
    elif pgaus == 2:
        x = [-0.577350269189625764509148780502, \
              0.577350269189625764509148780502]
        # ------------
        w = [1.0,\
             1.0]
    elif pgaus == 3:
        x = [-0.774596669241483377035853079956, \
              0.0, \
              0.774596669241483377035853079956]
        # ------------
        w = [5.0 / 9.0, \
             8.0 / 9.0, \
             5.0 / 9.0]
    elif pgaus == 4:
        x = [- 0.861136311594052575223946488893, \
             - 0.339981043584856264802665759103, \
               0.339981043584856264802665759103, \
               0.861136311594052575223946488893]
        # ------------
        w = [0.347854845137453857373063949222, \
             0.652145154862546142626936050778, \
             0.652145154862546142626936050778, \
             0.347854845137453857373063949222]
    elif pgaus == 5:
        x = [- 0.906179845938663992797626878299, \
             - 0.538469310105683091036314420700, \
               0.0, \
               0.538469310105683091036314420700, \
               0.906179845938663992797626878299]
        # -----------
        w = [0.236926885056189087514264040720, \
             0.478628670499366468041291514836, \
             0.568888888888888888888888888889, \
             0.478628670499366468041291514836, \
             0.236926885056189087514264040720]
    elif pgaus == 6:
        x = [- 0.932469514203152027812301554494, \
             - 0.661209386466264513661399595020, \
             - 0.238619186083196908630501721681, \
               0.238619186083196908630501721681, \
               0.661209386466264513661399595020, \
               0.932469514203152027812301554494]
        # -----------
        w = [0.171324492379170345040296142173, \
             0.360761573048138607569833513838, \
             0.467913934572691047389870343990, \
             0.467913934572691047389870343990, \
             0.360761573048138607569833513838, \
             0.171324492379170345040296142173]
    elif pgaus == 7:
        x = [- 0.949107912342758524526189684048, \
             - 0.741531185599394439863864773281, \
             - 0.405845151377397166906606412077, \
               0.0, \
               0.405845151377397166906606412077, \
               0.741531185599394439863864773281, \
               0.949107912342758524526189684048]
        # -----------
        w = [0.129484966168869693270611432679, \
             0.279705391489276667901467771424, \
             0.381830050505118944950369775489, \
             0.417959183673469387755102040816, \
             0.381830050505118944950369775489, \
             0.279705391489276667901467771424, \
             0.129484966168869693270611432679]
    elif pgaus == 8:
        x = [- 0.960289856497536231683560868569, \
             - 0.796666477413626739591553936476, \
             - 0.525532409916328985817739049189, \
             - 0.183434642495649804939476142360, \
               0.183434642495649804939476142360, \
               0.525532409916328985817739049189, \
               0.796666477413626739591553936476, \
               0.960289856497536231683560868569]
        # -----------
        w = [0.101228536290376259152531354310, \
             0.222381034453374470544355994426, \
             0.313706645877887287337962201987, \
             0.362683783378361982965150449277, \
             0.362683783378361982965150449277, \
             0.313706645877887287337962201987, \
             0.222381034453374470544355994426, \
             0.101228536290376259152531354310]
    elif pgaus == 9:
        x = [- 0.968160239507626089835576202904, \
             - 0.836031107326635794299429788070, \
             - 0.613371432700590397308702039341, \
             - 0.324253423403808929038538014643, \
               0.0, \
               0.324253423403808929038538014643, \
               0.613371432700590397308702039341, \
               0.836031107326635794299429788070, \
               0.968160239507626089835576202904]
        # -----------
        w = [0.812743883615744119718921581105E-01, \
             0.180648160694857404058472031243, \
             0.260610696402935462318742869419, \
             0.312347077040002840068630406584, \
             0.330239355001259763164525069287, \
             0.312347077040002840068630406584, \
             0.260610696402935462318742869419, \
             0.180648160694857404058472031243, \
             0.812743883615744119718921581105E-01]
    
    # Change type of arrays
    x = np.asarray(x); w = np.asarray(w)

    return x, w

def iga_find_positions_weights_element(degree, knotvector, el):
    " Computes Gauss weights and positions in isogeometric space for a known element "

    # Find knots of the knot-vector
    knots = np.unique(knotvector)

    # Find position and weight of Gauss points in isogeometric space
    xg, wg = gaussTable(degree + 1)
    
    # Find position of quadrature points at the element (scaling)
    xg = 0.5 * ((knots[el+1] - knots[el])*xg + knots[el] + knots[el+1])

    # Find weight of quadrature points at the element (scaling)
    wg = 0.5 * (knots[el+1] - knots[el]) * wg

    return xg, wg

def iga_find_positions_weights(degree, knotvector):
    " Computes Gauss weights and positions in parametric space for a known curve "

    # Find number of elements
    nbel = len(np.unique(knotvector)) - 1

    # Initialize vectors
    xg = []; wg = []

    for _ in range(nbel): 
        # Find quadrature pints for all elements
        xg_el, wg_el = iga_find_positions_weights_element(degree, knotvector, _)

        # Assmbly vectors
        xg.extend(xg_el); wg.extend(wg_el)

    return xg, wg

class IGA(thermoMechaModel):
    
    # ===========================
    # INITIALIZE 
    # ===========================

    def __init__(self, modelIGA: None, ElbyEl=False):
        super().__init__(modelIGA)

        if ElbyEl:
            # Set table of functions in every element (in each dimension)
            self.__set_table_funct_span()
        
            # Topology
            self.__set_bspline_topology()

            # Topology of quadrature points
            self.__set_qp_topology()

        # Quadrature points and weights in each dimension
        self.eval_positions_weights()

        # Eval basis in each dimension
        self.eval_basis()

        # Jacobien 
        self._Jqp, self._qp_PS = super().eval_jacobien_pps(self._dim, self._ctrlpts, self._DB, self._nb_qp_cgg_total)

        # Get conductivity and capacity coefficients
        self._conductivity_coef, self._capacity_coef, self._detJ = \
            super().eval_K_C_coefficient(self._nb_qp_cgg_total, self._Jqp, 
                                        self._conductivity, self._capacity)

        return

    def __set_table_funct_span(self):
        "Set table of functions on each element of the dimension"

        # Initiliaze
        self._table_funct_span = []

        for dim in range(self._dim):
            # Set table of functions per element 
            table = np.zeros((self._nb_el[dim][0], self._degree[dim][0] + 2), dtype= int); 
            table[0, 0] = self._degree[dim][0]; table[0, 1:] = np.arange(self._degree[dim][0] + 1) 

            # Set multiplicity
            multiplicity = 1

            for _ in range(1, self._nb_el[dim][0]): 
                # Set values of the table
                table[_, :2] = table[_-1, :2] + multiplicity
                table[_, 2:] = table[_, 1] + np.arange(1, self._degree[dim] + 1) 

            # Append
            self._table_funct_span.append(table)

        return
    
    def __set_bspline_topology(self): 
        """ Sets topology tables: 
        IEN: Element nodes
        INC: NURBS coordinates
        INE: NURBS elements
        """

        # Get number of dimensions
        dimensions = self._dim

        # Get number of control points in each dimension
        nb_ctrlpts = np.ones(3, dtype= int)
        for dim in range(dimensions):
            nb_ctrlpts[dim] = self._nb_ctrlpts[dim] 

        # Find total number of control points 
        nb_ctrlpts_total = self._nb_ctrlpts_total

        # Get number of control points in each dimension
        nb_el = np.ones(3, dtype= int)
        for dim in range(dimensions):
            nb_el[dim] = self._nb_el[dim] 

        # Find total number of control points 
        nb_el_total = self._nb_el_total

        # ----------------------
        # INC: NURBS coordinates
        # ----------------------
        INC = np.zeros((nb_ctrlpts_total, 3), dtype= int)

        for i3 in range(nb_ctrlpts[2]): 
            for i2 in range(nb_ctrlpts[1]): 
                for i1 in range(nb_ctrlpts[0]):
                    genPos = i1 + i2*(nb_ctrlpts[0]) + i3*(nb_ctrlpts[0]*nb_ctrlpts[1])
                    INC[genPos, :] = [i1, i2, i3] 

        # -------------------
        # INE: NURBS elements
        # -------------------
        INE = np.zeros((nb_el_total, 3), dtype= int)

        for i3 in range(nb_el[2]): 
            for i2 in range(nb_el[1]): 
                for i1 in range(nb_el[0]):
                    genPos = i1 + i2*(nb_el[0]) + i3*(nb_el[0]*nb_el[1])
                    INE[genPos, :] = [i1, i2, i3] 

        # -----------------------------
        # IEN: Element nodes
        # ----------------------------
        # Find IEN for each dimension
        IEN_dim = []
        for dim in range(self._dim): 
            IEN_dimt = np.zeros((self._nb_el[dim][0], self._degree[dim][0] + 1), dtype= int)
            IEN_dimt[0, :] = range(self._degree[dim][0] + 1)
        
            for _ in range(1, self._nb_el[dim][0]): 
                IEN_dimt[_, 0] = IEN_dimt[_-1, 0] + 1
                IEN_dimt[_, 1:] = IEN_dimt[_, 0] + np.arange(1, self._degree[dim][0] + 1) 

            IEN_dim.append(IEN_dimt.tolist())
        
        # Set global IEN
        IEN = []
        for el in range(self._nb_el_total):
            # Find position of global element in each dimension 
            el_dim = INE[el, :]

            # Find element nodes in each dimension
            EN_dim = []
            for dim in range(self._dim): 
                EN_dim.append(np.asarray(IEN_dim[dim])[el_dim[dim], :])
            
            # Find raw set of possible element nodes 
            set_dim = []    
            for dim in range(self._dim): 
                set_dimt = []
                for _ in EN_dim[dim]:
                    set_dimt.extend(np.where(INC[:, dim] == _)[0]) 
                set_dim.append(set_dimt)

            # Select element nodes
            EN = set(set_dim[0])
            for dim in range(1, self._dim): 
                EN = EN.intersection(set(set_dim[dim]))

            IEN.append(np.sort(list(EN)))

        # Assign values
        self._IEN, self._INC, self._INE = IEN, INC, INE

        return

    def __set_qp_topology(self):
        """ Sets topology tables: 
        IEN: Table of quadrature points nodes
        """

        # Get number of dimensions
        dimensions = self._dim

        # Get number of control points in each dimension
        nb_qp_cgg = np.ones(3, dtype= int)
        for dim in range(dimensions):
            nb_qp_cgg[dim] = self._nb_qp_cgg[dim][0]

        # Find total number of control points 
        nb_qp_cgg_total = self._nb_qp_cgg_total
        
        # ---------------------------------
        # INC adapted for quadrature points
        # ---------------------------------
        # Initialize table
        INC = np.zeros((nb_qp_cgg_total, 3), dtype= int)

        for i3 in range(nb_qp_cgg[2]): 
            for i2 in range(nb_qp_cgg[1]): 
                for i1 in range(nb_qp_cgg[0]):
                    genPos = i1 + i2*(nb_qp_cgg[0]) + i3*(nb_qp_cgg[0]*nb_qp_cgg[1])
                    INC[genPos, :] = [i1, i2, i3] 

        # ---------------------------------
        # INE adapted for quadrature points
        # ---------------------------------
        INE = self._INE

        # -----------------------------
        # IEN adapted for quadrature points
        # ----------------------------
        # Find IEN for each dimension
        IEN_dim = []
        for dim in range(self._dim): 
            IEN_dimt = np.zeros((self._nb_el[dim][0], self._degree[dim][0] + 1), dtype= int)
            IEN_dimt[0, :] = range(self._degree[dim][0] + 1)
        
            for _ in range(1, self._nb_el[dim][0]): 
                IEN_dimt[_, 0] = IEN_dimt[_-1, -1] + 1
                IEN_dimt[_, 1:] = IEN_dimt[_, 0] + np.arange(1, self._degree[dim][0] + 1) 

            IEN_dim.append(IEN_dimt.tolist())
        
        # Set global IEN adapted for quadrature points
        IEN = []

        for el in range(self._nb_el_total):
            # Find position of global element in each dimension 
            el_dim = INE[el, :]

            # Find quadrature points nodes in each dimension
            QPN_dim = []
            for dim in range(self._dim): 
                QPN_dim.append(np.asarray(IEN_dim[dim])[el_dim[dim], :])
            
            # Find raw set of possible quadrature points nodes 
            set_dim = []    
            for dim in range(self._dim): 
                set_dimt = []
                for _ in QPN_dim[dim]:
                    set_dimt.extend(np.where(np.asarray(INC[:, dim]) == _)[0]) 
                set_dim.append(set_dimt)

            # Select quadrature points nodes
            QPN = set(set_dim[0])
            for dim in range(1, self._dim): 
                QPN = QPN.intersection(set(set_dim[dim]))

            IEN.append(np.sort(list(QPN)))

        # Assign values
        self._IEN_qp, self._INC_qp = IEN, INC

        return IEN, INC

    def eval_positions_weights(self):
        " Computes weights and position of quadrature points "

        # Initialize quadrature weigths and positions
        self._qp_cgg_dim = []
        self._DW = []

        for dim in range(self._dim): 
            qp_dim_temp, weight_dim_temp = iga_find_positions_weights(self._degree[dim][0], self._knotvector[0][dim])
            self._qp_cgg_dim.append(qp_dim_temp)
            self._DW.append(weight_dim_temp)

        return
    
    def eval_basis(self):
        " Computes basis at quadrature points "

        print('Evaluating basis')
        start = time.time()
        self._DB = []
        for dim in range(self._dim): 
            B0, B1 = eval_basis(self._degree[dim][0], self._knotvector[0][dim], self._qp_cgg_dim[dim])
            self._DB.append([B0, B1]) 
        stop = time.time()
        print('\tBasis in : %.5f s' %(stop-start))
        
        return

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        source_coef = super().eval_F_coefficient(fun, self._dim, self._qp_PS, self._detJ)
        return source_coef

    # ===========================
    # ASSEMBLE 
    # ===========================

    def eval_conductivity_matrix(self):
        " Assemble conductivity matrix K "

        start = time.time()

        # Initialize conductivity matrix 
        K = sp.csr_matrix((self._nb_ctrlpts_total, self._nb_ctrlpts_total))

        Wt = 1
        for dim in range(self._dim):
            Wt = sp.kron(sp.diags(self._DW[dim]), Wt)

        for j in range(self._dim):
            beta = np.zeros(self._dim, dtype = int); beta[j] = 1
            Btr = 1 
            for dim in range(self._dim):
                bt = beta[dim]
                Btr = sp.kron(self._DB[dim][bt], Btr)
            
            for i in range(self._dim):
                Cij = self._conductivity_coef[i, j, :]
                alpha = np.zeros(self._dim, dtype = int); alpha[i] = 1
                
                Btl = 1
                for dim in range(self._dim):
                    at = alpha[dim] 
                    Btl = sp.kron(self._DB[dim][at], Btl)
                    
                # Evaluates Cij * B in each point
                Btr_Cij = sp.csr_matrix.dot(Btr, sp.diags(Cij))

                # Find K
                BW = sp.csr_matrix.dot(Btl.tocsr()[:,:], Wt.tocsr()[:,:])
                K += sp.csr_matrix.dot(BW.tocsr()[:,:], Btr_Cij.tocsr()[:,:].T)

        stop = time.time()
        print('Conductivity matrix assembled in : %.5f s' %(stop-start))

        return K

    def eval_capacity_matrix(self):
        " Assemble capacity matrix C "

        start = time.time()
        
        # Initialize weights
        B = 1; W = 1
        for dim in range(self._dim): 
            # Find basis
            B0 = self._DB[dim][0]

            # Find weights
            W00 = self._DW[dim]

            # Find W and B
            B = sp.kron(B0, B)
            W = np.kron(np.array(W00), W)
            
        # Assemble C
        W = sp.csr_matrix.dot(W, sp.diags(self._capacity_coef))
        C = sp.csr_matrix.dot(sp.csr_matrix.dot(B.tocsr()[:,:], sp.diags(W)), B.tocsr()[:,:].transpose())

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))

        return C

    def eval_source_vector(self, fun):
        " Assemble power density vector F "

        start = time.time()
        
        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)

        # Initialize weights
        B = 1; W = 1
        for dim in range(self._dim): 
            # Find basis
            B0 = self._DB[dim][0]

            # Find weights
            W00 = self._DW[dim]

            # Find W and B
            W = sp.kron(sp.diags(W00), W)
            B = sp.kron(B0, B)

        # Assemble F
        F = sp.csr_matrix.dot(sp.csr_matrix.dot(B.tocsr()[:,:], W.tocsr()[:,:]), self._source_coef)

        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return F

    def eval_conductivity_matrix_element(self):
        " Assemble conductivity matrix K element by element "

        start = time.time()

        # Initialize conductivity matrix 
        K = sp.lil_matrix((self._nb_ctrlpts_total, self._nb_ctrlpts_total))

        for el in range(self._nb_el_total): 
            # Find NURBS elements of global element 
            INE_el = self._INE[el, :]

            # Intialize basis and weights
            DBel = []; DWel = []

            for dim in range(self._dim):
                # Find NURBS element position in each dimension 
                el_dim = INE_el[dim]

                # Set indexes (positions) of quadrature points in each dimension  
                ind_xg = np.arange(el_dim*(self._degree[dim][0]+1), (el_dim+1)*(self._degree[dim][0]+1), dtype= int)

                # Find indexes of functions in NURBS element in each dimension
                ind_basis = self._table_funct_span[dim][el_dim][1:]

                # Select weights and basis in the element 
                wgel = np.asarray(self._DW[dim])[ind_xg]
                B0el = self._DB[dim][0].tocsr()[np.ix_(ind_basis, ind_xg)]
                B1el = self._DB[dim][1].tocsr()[np.ix_(ind_basis, ind_xg)]

                # Save 
                DWel.append(wgel)
                DBel.append([B0el, B1el])

            # Find element nodes in global element
            IEN_el = np.asarray(self._IEN)[el, :].tolist()

            # Find quadrature points in element
            IEN_QP_el = np.asarray(self._IEN_qp)[el, :]

            # Select multiplier of C
            conductivity_coef_el = np.zeros((self._dim, self._dim, len(IEN_QP_el)))
            conductivity_coef_el[:, :, :] = self._conductivity_coef[:, :, IEN_QP_el]
            
            # Initialize elementary conductivity matrix 
            Kel = sp.lil_matrix((len(IEN_el), len(IEN_el)))

            for i in range(self._dim):
                for j in range(self._dim):
                    Cij = conductivity_coef_el[i, j, :]
                    
                    alpha = np.zeros(self._dim, dtype = int); alpha[i] = 1
                    beta = np.zeros(self._dim, dtype = int); beta[j] = 1

                    Btl = 1; Btr = 1; Wt = 1
                    for dim in range(self._dim):
                        at = alpha[dim] 
                        bt = beta[dim]

                        Btl = sp.kron(DBel[dim][at], Btl)
                        Btr = sp.kron(DBel[dim][bt], Btr)
                        Wt = sp.kron(sp.diags(DWel[dim]), Wt)

                    # Evaluates Cij * B in each point
                    Btr = sp.csr_matrix.dot(Btr, sp.diags(Cij))

                    # Find elementary K
                    BW = sp.csr_matrix.dot(Btl.tocsr()[:,:], Wt.tocsr()[:,:])
                    Kel += sp.csr_matrix.dot(BW.tocsr()[:,:], Btr.tocsr()[:,:].T)

            # Assemble of K
            K[np.ix_(IEN_el, IEN_el)] += Kel

        stop = time.time()
        print('Conductivity matrix assembled in : %.5f s' %(stop-start))

        return K

    def eval_capacity_matrix_element(self):
        " Assemble capacity matrix C element by element "

        start = time.time()

        # Initialize capacity matrix 
        C = sp.lil_matrix((self._nb_ctrlpts_total, self._nb_ctrlpts_total))

        for el in range(self._nb_el_total): 
            # Find NURBS elements of global element 
            INE_el = self._INE[el, :]

            # Intialize basis and weights
            B = 1; W = 1

            for dim in range(self._dim): 
                # Find NURBS element position in each dimension 
                el_dim = INE_el[dim]

                # Set indexes (positions) of quadrature points in each dimension  
                ind_xg = np.arange(el_dim*(self._degree[dim][0]+1), (el_dim+1)*(self._degree[dim][0]+1), dtype= int)

                # Find indexes of functions in NURBS element in each dimension
                ind_basis = self._table_funct_span[dim][el_dim][1:]

                # Select weights and basis in the element 
                wgel = np.asarray(self._DW[dim])[ind_xg]
                B0el = self._DB[dim][0].tocsr()[np.ix_(ind_basis, ind_xg)]
                
                # Find W and B
                W = sp.kron(sp.diags(wgel), W)
                B = sp.kron(B0el, B)

            # Find element nodes in global element
            IEN_el = np.asarray(self._IEN)[el, :].tolist()

            # Find quadrature points in element
            IEN_QP_el = np.asarray(self._IEN_qp)[el]

            # Select multiplier of C
            capacity_coef_el = self._capacity_coef[IEN_QP_el] 

            # Evaluate B W
            BW = sp.csr_matrix.dot(B.tocsr()[:,:], W.tocsr()[:,:])

            # Evaluate B
            B = B * sp.diags(capacity_coef_el)
            
            # Aseemble C
            Cel = sp.csr_matrix.dot(BW.tocsr()[:,:], B.tocsr()[:,:].T)

            C[np.ix_(IEN_el, IEN_el)] += Cel

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))

        return C

    def eval_source_vector_element(self, fun): 
        " Assemble power density vector F element by element "

        start = time.time()

        # Get dimension
        dim = self._dim
      
        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)

        # Initialize heat vector 
        F = np.zeros(self._nb_ctrlpts_total)

        for el in range(self._nb_el_total): 
            # Find NURBS elements of global element 
            INE_el = self._INE[el, :]

            # Intialize basis and weights
            B = 1; W = 1

            for dim in range(self._dim): 
                # Find NURBS element position in each dimension 
                el_dim = INE_el[dim]

                # Set indexes (positions) of quadrature points in each dimension  
                ind_xg = np.arange(el_dim*(self._degree[dim][0]+1), (el_dim+1)*(self._degree[dim][0]+1), dtype= int)

                # Find indexes of functions in NURBS element in each dimension
                ind_basis = self._table_funct_span[dim][el_dim][1:]

                # Select 
                wgel = np.asarray(self._DW[dim])[ind_xg]
                B0el = self._DB[dim][0].tocsr()[np.ix_(ind_basis, ind_xg)]

                # Find W and B
                W = sp.kron(sp.diags(wgel), W)
                B = sp.kron(B0el, B)

            # Find element nodes in global element
            IEN_el = np.asarray(self._IEN)[el, :].tolist()

            # Find quadrature points in element
            IEN_QP_el = np.asarray(self._IEN_qp)[el]

            # Find elementary F
            source_coef_el = np.asarray(self._source_coef)[IEN_QP_el] 

            # Assembly of F
            Fel = sp.csr_matrix.dot(sp.csr_matrix.dot(B.tocsr()[:,:], W.tocsr()[:,:]), source_coef_el)
            F[IEN_el] += Fel

        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return F
