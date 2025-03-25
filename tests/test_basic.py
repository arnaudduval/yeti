def test_imports():
    import yeti_iga as yeti

    from yeti_iga.preprocessing.igaparametrization import IGAparametrization
    from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip
    from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
    import yeti_iga.reconstructionSOL as rsol
    import yeti_iga.postprocessing.postproc as pp