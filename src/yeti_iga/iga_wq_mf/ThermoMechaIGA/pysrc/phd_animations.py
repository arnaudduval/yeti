"""
.. Isogeometric analysis animations 
.. Joaquin Cornejo
.. How to use: in terminal call "manim -pql filename classname" 
"""

# Python libraries
from manim import *
from scipy import sparse
import cv2

# My libraries
from lib.base_functions import create_knotvector, eval_basis_python, wq_find_positions

class BSplineCurve(Scene):
    def construct(self):
        # Define color
        colors = [BLUE, TEAL, GREEN, YELLOW, GOLD, RED, MAROON, PURPLE, PINK]

        # Create first graph
        ax1 = Axes(
            x_range=[0, 1, 0.25],
            y_range=[-0.1, 1.1, 0.2],
            x_length=4,
            y_length=3,
            x_axis_config={"numbers_to_include": np.arange(0, 1.1, 0.25)},
            y_axis_config={"numbers_to_include": np.arange(0, 1.1, 0.2)},
            tips=False,
        )
        ax1.shift(LEFT*2.5)

        # Set labels of first graph
        labels1 = ax1.get_axis_labels(x_label=Tex(r"$\xi$"), y_label=Tex(r"$B_i(\xi)$"))
        title1 = Text("Parametric space", font_size= 20).move_to(LEFT*2.5+DOWN*2.1)
        self.add(ax1, labels1, title1)

        # Create second graph
        ax2 = Axes(
            x_range=[-1.1, 2.1, 0.5],
            y_range=[-1.1, 2.1, 0.5],
            x_length=4,
            y_length=3,
            x_axis_config={"numbers_to_include": np.arange(-1, 2, 1)},
            y_axis_config={"numbers_to_include": np.arange(-1, 2, 1)},
            tips=False,
        )
        ax2.shift(RIGHT*3)
        # Set labels of first graph
        labels2 = ax2.get_axis_labels(x_label=Tex("x"), y_label=Tex("y"))
        title2 = Text("Physical space", font_size= 20).move_to(RIGHT*3+DOWN*2.1)

        # Define B-Spline functions
        # -------------------------
        N = 100
        degree = 2
        nbel = 4
        knots = np.linspace(0, 1, N)
        knotvector = create_knotvector(degree, nbel)
        B0, ifunt = eval_basis_python(degree, knotvector, knots)
        B0 = B0.toarray()

        # Define B-spline curve
        # ---------------------
        ctrlpts = np.asarray([[-1, 1, 0], [-0.5, 0.25, 0], [0, 2, 0], 
                    [0.75, -0.5, 0], [1.5, 1, 0], [2, 0, 0]])

        graphDotSet = ax2.plot_line_graph(x_values=ctrlpts[:,0], y_values=ctrlpts[:,1], 
                                        add_vertex_dots=True, line_color= BLUE_A, stroke_width = 2)
        self.add(ax2, labels2, title2, graphDotSet)

        # Draw curve
        xx = np.zeros(N)
        yy = np.zeros(N)
        for i, ifunt in enumerate([0, 5, 1, 4, 2, 3]):
            # Plot Bezier functions
            B0_knots = B0[ifunt, :]
            graph1 = ax1.plot_line_graph(x_values=knots, y_values=B0_knots, 
                                        add_vertex_dots=False, line_color= colors[ifunt])

            # Plot B-spline curve
            xx += ctrlpts[ifunt][0]*B0[ifunt, :]
            yy += ctrlpts[ifunt][1]*B0[ifunt, :]
            graph2 = ax2.plot_line_graph(x_values=xx, y_values=yy, 
                                        add_vertex_dots=False, line_color= PINK)

            graphDot = Dot(ax2.coords_to_point(ctrlpts[ifunt][0], ctrlpts[ifunt][1]), color=colors[ifunt])
            
            # Animate
            self.add(graphDot)
            self.play(Create(graph1), Create(graph2))    
            self.wait(1)
            if i != degree+nbel-1: self.play(FadeOut(graph2))
            del graph1, graph2

        return

class WQExplication(Scene):
    def construct(self):
        # Define color
        colors = [BLUE, TEAL, GREEN, YELLOW, GOLD, RED, MAROON, PURPLE, PINK]

        # Create first graph 
        ax0 = Axes(
            x_range=[0, 1, 0.25],
            y_range=[-0.1, 0.1, 0.2],
            x_length=4,
            y_length=2,
            x_axis_config={"numbers_to_include": np.arange(0, 1.1, 0.25)},
            tips=False,
        )
        # Set labels of first graph
        labels0 = ax0.get_axis_labels(x_label=Tex(r"$\xi$"), y_label=Tex(""))
        groupe0 = VGroup(ax0, labels0)

        # Create second graph
        ax1 = Axes(
            x_range=[0, 1, 0.25],
            y_range=[-0.1, 1.1, 0.2],
            x_length=4,
            y_length=3,
            x_axis_config={"numbers_to_include": np.arange(0, 1.1, 0.25)},
            y_axis_config={"numbers_to_include": np.arange(0, 1.1, 0.2)},
            tips=False,
        )
        # Set labels of second graph
        labels1 = ax1.get_axis_labels(x_label=Tex(r"$\xi$"), y_label=Tex(r"$B_i(\xi)$"))
        groupe1 = VGroup(ax1, labels1)

        # Define B-spline functions
        # -------------------------
        N = 100
        degree = 2
        nbel = 4
        knots = np.linspace(0, 1, N)
        knotvector = create_knotvector(degree, nbel)
        B0, _ = eval_basis_python(degree, knotvector, knots)
        B0 = B0.toarray()

        # Define quadrature points
        qp_wq = wq_find_positions(degree, knotvector, 2)
        B0_wq, _ = eval_basis_python(degree, knotvector, qp_wq)
        B0_wq = B0_wq.toarray()

        # Scene 1
        # ----------------------
        groupe0.add(ax0.plot_line_graph(x_values=qp_wq, y_values=np.zeros(len(qp_wq)), 
                                        add_vertex_dots=True, line_color= WHITE, stroke_width = 2))
        groupe1.add(ax1.plot_line_graph(x_values=qp_wq, y_values=np.zeros(len(qp_wq)), 
                                        add_vertex_dots=True, line_color= WHITE, stroke_width = 2))
        self.play(Create(groupe0))
        self.wait(1)
        self.play(ReplacementTransform(groupe0, groupe1))
        for _ in range(degree+nbel):
            # Plot Bezier functions
            B0_knots = B0[_, :]
            graph1 = ax1.plot_line_graph(x_values=knots, y_values=B0_knots, 
                                        add_vertex_dots=False, line_color= colors[_])
            groupe1.add(graph1)
            self.play(Create(graph1), run_time=0.2)  
        self.wait(1)
        self.play(groupe1.animate.shift(3*LEFT).scale(0.8))
        self.wait(1)
         
        # Scene 2
        # ----------------------
        # Define matrices and vectors
        vm0 = Matrix([["B_1"], ["B_2"], ["B_3"]])         
        m0 = Matrix([[1., 0.44,  0.11], 
                    [0., 0.5, 0.67], 
                    [0., 0.055,	0.22]], h_buff=1.5)
        equalsymbol0 = Tex("=")

        vm1 = Matrix([["I_{11}"], ["I_{12}"], ["I_{13}"]], 
            left_bracket="\{",
            right_bracket="\}")         
        m1 = Matrix([[0.05], [0.0292], [0.0042]], 
            left_bracket="\{",
            right_bracket="\}")
        equalsymbol1 = Tex("=")

        self.play(
            Write(vm0.move_to(RIGHT*0.25+UP*2.5)),
            Write(equalsymbol0.move_to(RIGHT*1.25+UP*2.5)),
            Write(m0.move_to(RIGHT*4+UP*2.5)))
        self.wait()

        self.play(
            Write(vm1.move_to(RIGHT*0.25+DOWN*1)),
            Write(equalsymbol1.move_to(RIGHT*1.75+DOWN*1)),
            Write(m1.move_to(RIGHT*4+DOWN*1)))
        self.wait()

        equation = MathTex(
            "\\mathbb{B}","\\mathbb{W}","=",
            "\\mathbb{I}"
        )
        self.play(Write(equation.move_to(RIGHT*2+DOWN*3)))

        # For equation
        framebox1 = SurroundingRectangle(equation[0], buff = .1)
        framebox2 = SurroundingRectangle(equation[3], buff = .1)

        # For matrix
        framebox3 = SurroundingRectangle(m0.get_rows(), buff = .1)
        framebox4 = SurroundingRectangle(m1.get_rows(), buff = .1)

        self.play(
            Create(framebox1),
            Create(framebox3),
        )
        self.wait()
        self.play(
            ReplacementTransform(framebox1,framebox2),
            ReplacementTransform(framebox3,framebox4),
        )
        self.wait(2)

        return

class BSplineSurface(Scene):
    def construct(self):
        # Define color
        colors = [BLUE, TEAL, GREEN, YELLOW, GOLD, MAROON, PURPLE, PINK]

        # Create first graph
        ax1 = Axes(
            x_range=[0, 1, 0.25],
            y_range=[0, 1, 0.25],
            x_length=3,
            y_length=3,
            x_axis_config={"numbers_to_include": np.arange(0, 1.1, 0.25)},
            y_axis_config={"numbers_to_include": np.arange(0, 1.1, 0.25)},
            tips=False,
        )
        ax1.shift(RIGHT*1+DOWN*1)
        # Set labels of first graph
        labels1 = ax1.get_axis_labels(x_label=Tex(r"$\xi$"), y_label=Tex(r"$\nu$"))

        # Create second graph
        ax2 = Axes(
            x_range=[0, 1, 0.25],
            y_range=[0, 1, 0.25],
            x_length=3,
            y_length=1.5,
            axis_config={"include_numbers": False},
            tips=False,
        )
        ax2.shift(RIGHT*1+UP*3)
        # Set labels of first graph
        labels2 = ax2.get_axis_labels(x_label=Tex(r"$\xi$"), y_label=Tex(r"$B_i(\xi)$"))

        # Create third graph
        ax3 = Axes(
            x_range=[0, 1, 0.25],
            y_range=[0, 1, 0.25],
            x_length=3,
            y_length=1.5,
            axis_config={"include_numbers": False},
            tips=False,
        )
        ax3.shift(LEFT*3+DOWN*1)
        # Set labels of first graph
        labels3 = ax3.get_axis_labels(x_label=Tex(r"$\nu$"), y_label=Tex(r"$B_i(\nu)$"))
        ax3.rotate(PI/2)
        labels3.rotate(PI/2)

        # Define Bezier surface
        # -------------------
        N = 256
        degree = 3
        nbel = 1
        knots = np.linspace(0, 1, N)
        knotvector = create_knotvector(degree, nbel)
        B0, _ = eval_basis_python(degree, knotvector, knots)
        B0 = B0.toarray()
        B = sparse.kron(B0, B0)

        self.add(ax1, labels1)
        self.add(ax2, labels2)
        self.add(ax3, labels3)

        for n in range(degree+nbel):
            B_nu = B0[n, :]
            graph2 = ax3.plot_line_graph(x_values=knots, y_values=B_nu, add_vertex_dots=False, line_color= colors[n])
            self.add(graph2)

            for m in range(degree+nbel):
                B_xi = B0[m, :]
                graph1 = ax2.plot_line_graph(x_values=knots, y_values=B_xi, add_vertex_dots=False, line_color= colors[m])
                self.add(graph1)

                # Control point using
                ctrlpts_pos = m + n * (degree+nbel)
                B0_2D = B.getrow(ctrlpts_pos).toarray().reshape((N, N))
                B0_2D = np.flipud(np.uint8(B0_2D*255))
                im_color = cv2.applyColorMap(B0_2D, cv2.COLORMAP_INFERNO)

                # Objects
                image_par = ImageMobject(im_color, image_mode="RBG").scale(1.57)
                image_par.set_resampling_algorithm(RESAMPLING_ALGORITHMS["cubic"])
                image_par.shift(RIGHT*1+DOWN*1)

                self.play(FadeIn(image_par), FadeIn(graph1), run_time=0.2)
                self.wait(0.5)
                self.play(FadeOut(image_par), FadeOut(graph1), run_time=0.2)

            self.play(FadeOut(graph2), run_time=0.5)

        return

class SumFactorization(Scene):
    def construct(self):
        # Scene 1
        # ----------------------
        # Define equation
        equation1 = MathTex(
            "(","V^{(2)}","\\otimes","V^{(1)}", ")", "x")
        equation1[1].set_color(GREEN)
        equation1[3].set_color(BLUE)
        equation1[5].set_color(MAROON)
        self.play(Write(equation1.move_to(UP*3)))
        text = Tex("Using Kronecker product :")
        self.play(Write(text.move_to(UP*2+LEFT*2)))

        # Defines vectors and matrices    
        v1 = Matrix([[1, 0, 2]], h_buff=0.6).set_color(BLUE)
        v1.get_brackets().set_color(BLUE).scale(0.7)

        v2 = Matrix([[2, 3]], h_buff=0.6).set_color(GREEN)
        v2.get_brackets().set_color(GREEN).scale(0.7)

        vr1 = Matrix([[2, 0, 4, 3, 0, 6]], v_buff=0.5).scale(0.7)
        vr2 = Matrix([[12]], v_buff=0.5).scale(0.7)

        x = Matrix([[1], [2], [1], [0], [2], [1]],
                    v_buff=0.7).set_color(MAROON).scale(0.7)
        x.get_brackets().set_color(MAROON)

        Xtensor = Matrix([[1, 0], [2, 2], [1, 1]],
                    v_buff=0.7).set_color(MAROON).scale(0.7)
        Xtensor.get_brackets().set_color(MAROON)

        # Define symbols
        kronprod = MathTex("\\otimes")
        dotprod0 = MathTex("\\dot")
        dotprod1 = MathTex("\\dot")

        g1 = VGroup(v2, kronprod, v1, dotprod0, x).arrange(buff=1)
        g2 = VGroup(vr1, dotprod1, x).arrange(buff=1)
        self.play(Create(g1))
        self.wait()
        self.play(ReplacementTransform(g1, g2))
        self.wait()
        self.play(ReplacementTransform(g2, vr2))
        self.wait()
        self.play(FadeOut(vr2), FadeOut(text))
        del  g1, g2, v1, v2, x
        # ============================

        # Scene 2
        # ----------------------
        text = Tex("Using sum-factorization algorithm :")
        self.play(Write(text.move_to(UP*2+LEFT*2)))

        # Define vectors and matrices     
        v1 = Matrix([[1, 0, 2]], h_buff=0.6).set_color(BLUE)
        v1.get_brackets().set_color(BLUE).scale(0.7)

        v2 = Matrix([[2, 3]], h_buff=0.6).set_color(GREEN)
        v2.get_brackets().set_color(GREEN).scale(0.7)

        x = Matrix([[1], [2], [1], [0], [2], [1]],
                    v_buff=0.7).set_color(MAROON).scale(0.7)
        x.get_brackets().set_color(MAROON)

        g1 = VGroup(v2, v1, x).arrange(buff=1)
        g2 = VGroup(v2, v1, Xtensor).arrange(buff=1)
       
        self.play(Create(g1))
        self.wait()

        #  Define equation
        equation2 = MathTex(
            "(","V^{(2)}","\\otimes","V^{(1)}", ")", 
            "x", "=", "\\sum_{j_1=1}", "\\sum_{j_2=1}", "V^{(2)}_{j_2}", "V^{(1)}_{j_1}", 
            "\\mathcal{X}_{j_1 j_2}"
        )
        equation2[1].set_color(GREEN)
        equation2[3].set_color(BLUE)
        equation2[5].set_color(MAROON)
        equation2[9].set_color(GREEN)
        equation2[10].set_color(BLUE)
        equation2[11].set_color(MAROON)
        self.play(ReplacementTransform(equation1, equation2.move_to(UP*3)))
        self.wait()
        self.play(ReplacementTransform(g1, g2))
        self.play(g2.animate.move_to(UP*0.5))
        self.wait()

        # Scene 3
        # ----------------------
        v1rectangle_old = SurroundingRectangle(v1.get_columns()[0])
        v2rectangle_old = SurroundingRectangle(v2.get_columns()[0])
        Xrectangle_old = SurroundingRectangle(Xtensor.get_entries()[0])
        text_algo_old = Tex("For $j_2$ = " + str(1) + " and " + "$j_1$ = " + str(1) + " :")
        self.play(Create(text_algo_old.move_to(DOWN)))

        # The first iteration
        r0, r1, result = 0, 0, 0
        result_text_old = Tex("S = " + str(r0) + " + " + str(r1) + " = " + str(result))
        self.play(Create(result_text_old.move_to(DOWN*2)))
        
        # Example values
        v1ex = [1, 0, 2]
        v2ex = [2, 3]
        xex = [1, 2, 1, 0, 2, 1]
        Xtex = np.array(xex).reshape((len(v1ex), len(v2ex)), order='F')
        for j2 in range(2):
            for j1 in range(3):
                # Update new values
                r1 = v1ex[j1] * v2ex[j2] * Xtex[j1, j2] 
                result += r1

                text_algo = Tex("For $j_2$ = " + str(j2+1) + " and " + "$j_1$ = " + str(j1+1) + " :")
                v1rectangle = SurroundingRectangle(v1.get_columns()[j1])
                v2rectangle = SurroundingRectangle(v2.get_columns()[j2])
                Xrectangle = SurroundingRectangle(Xtensor.get_entries()[j2+j1*2])
                result_text = Tex("S = " + str(r0) + " + " + str(r1) + " = " + str(result))

                # Animations
                self.play(
                    ReplacementTransform(v1rectangle_old,v1rectangle),
                    ReplacementTransform(v2rectangle_old,v2rectangle), 
                    ReplacementTransform(Xrectangle_old,Xrectangle), 
                    ReplacementTransform(text_algo_old.move_to(DOWN),text_algo.move_to(DOWN)),
                    ReplacementTransform(result_text_old.move_to(DOWN*2),result_text.move_to(DOWN*2))
                )
                
                # Update old values
                v1rectangle_old = v1rectangle
                v2rectangle_old = v2rectangle
                Xrectangle_old = Xrectangle
                text_algo_old = text_algo
                result_text_old = result_text
                r0 = result
                
                self.wait()

        return

class FastDiag(Scene):
    def construct(self):
        
        # Create sector
        sector = Sector(outer_radius=2, inner_radius=1).move_to(4*LEFT)
        sector.set_color(ORANGE).scale(1.2)
        self.play(Create(sector))

        text1 = MarkupText("Heat transfer equation is:", 
                    font_size=40, color=ORANGE).move_to(RIGHT*1+UP*1)
        equation1 = MathTex(
            "\\mathbb{K} \cdot T =  \\mathbb{F}").move_to(RIGHT*1)
        text2 = Tex("Stiffness matrix "+r"$\mathbb{K}$"+"\\\\includes geometric and thermal properties.", 
                    font_size=35).move_to(RIGHT*1+DOWN*1)
        self.play(Write(text1), Write(equation1), Write(text2))
        self.wait(2)
        text3 = Tex("To improve condition number of " + r"$\mathbb{K}$"+",\\\\we use a preconditioner " + r"$\mathcal{P}$",
                    font_size=35).move_to(RIGHT*1+DOWN*1)
        equation2 = MathTex(
            "\\mathcal{P}^{-1} \\cdot \\mathbb{K} \cdot T =  \\mathcal{P}^{-1} \\cdot \\mathbb{F}").move_to(RIGHT*1)
        self.play(ReplacementTransform(text2, text3))
        self.wait(2)
        self.play(ReplacementTransform(equation1, equation2))    
        self.wait(2)
        self.play(FadeOut(text3, text1, equation2))
        del text1, text3, equation2, equation1, text2

        # ============================================
        square = Square(side_length=2, color=BLUE).move_to(4*LEFT)
        square.set_fill(BLUE, opacity=0.5) .scale(1.2)
        self.play(ReplacementTransform(sector, square))

        text1 = MarkupText("Fast diagonalization:", 
                    font_size=40, color=BLUE).move_to(RIGHT*2+UP*1)
        text2 = Tex(r"$\mathcal{P}$"+" is the stiffness matrix of a square\\\\with thermal properties equal to the identity", 
                    font_size=35).move_to(RIGHT*2)
        self.play(Write(text1), Write(text2))
        self.wait(1)
        equation1 = MathTex(
            "\\mathcal{P}_{i, j} = \\int_{\\Omega} \\nabla B_i(\\hat{x}) \, \\nabla B_j(\\hat{x})\, d\\hat{x}", 
            ).move_to(RIGHT*1.75+DOWN*1)
        self.play(Write(equation1))
        self.wait(2)

        equation2 = MathTex(
            "\\mathcal{P} = K_2 \\otimes M_1 + M_2 \\otimes K_1 \\\\"
        ).move_to(RIGHT*1.75+DOWN*1)

        equation3 = MathTex(
            "K = \\int_0^1 B'(\\xi) B'(\\xi) d\\xi \\quad M = \\int_0^1 B(\\xi) B(\\xi) d\\xi"
        ).scale(0.7).move_to(RIGHT*1.75+DOWN*2)

        self.play(ReplacementTransform(equation1, equation2), Write(equation3))
        self.wait(2)
        self.play(FadeOut(equation3, text2), ApplyMethod(equation2.shift,1*UP))
        del equation3, text2, equation1
        self.wait(2)

        # =======================================
        text2 = Tex("Applying eigendecomposition \\\\and tensor product properties:", font_size=35).move_to(RIGHT*2)
        equation1 = MathTex(
            "\\mathcal{P} = (U_2\\otimes U_1)^{-T} (I_2 \\otimes D_1 + D_2 \\otimes I_1) (U_2\\otimes U_1)^{-1}"
        ).move_to(RIGHT*2+DOWN*1).scale(0.75)
        text3 = Tex("Computing " + r"$s = \mathcal{P}^{-1} r$" + " becomes feasible", font_size=35).move_to(RIGHT*2+DOWN*2)
        self.play(Write(text2), ReplacementTransform(equation2, equation1), Write(text3))
        self.wait(3)

        return