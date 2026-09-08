class UpdateManager:
    def __init__(self):
        self.step_count = 0
        self.newton_iter = 0

    def reset_newton(self):
        self.newton_iter = 0

    def increment_newton(self):
        self.newton_iter += 1

    def increment_step(self):
        self.step_count += 1
        self.reset_newton()

    @property
    def should_update_material(self):
        # Rule: always update (every newton in every time step)
        # Eventually, we can adapt more complicated rules
        return True

    @property
    def should_update_preconditioner(self):
        # Rule: Update at the start of a time step
        # Eventually, we can adapt more complicated rules
        return self.newton_iter == 0

    @property
    def should_update_penalty(self):
        # Rule: never update penalty (except first newton at first time step)
        # Or when preconditioner should be updated
        # That is because Lagrange Preconditioner will also depends
        # on Fast Diagonalization preconditioner
        # Eventually, we can adapt more complicated rules
        return (
            self.step_count == 0 and self.newton_iter == 0
        ) or self.should_update_preconditioner

    def get_status(self):
        return f"UpdateManager, step: {self.step_count}, iteration: {self.newton_iter}"
