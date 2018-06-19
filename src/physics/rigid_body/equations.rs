use integrate::RK2;
use super::RigidBody;

impl RK2 for RigidBody {
    fn initialize(&mut self, dt: f32){

    }
    fn stage1(&mut self, dt: f32){
        let dtb2 = dt / 2.;
        self.r_com[0] = self.r_com[0] + dtb2 * self.v_com[0];
        self.r_com[1] = self.r_com[1] + dtb2 * self.v_com[1];
        self.r_com[2] = self.r_com[2] + dtb2 * self.v_com[2];

        self.v_com[0] = self.v_com[0] + dtb2 * self.frc[0] / self.m_total;
        self.v_com[1] = self.v_com[1] + dtb2 * self.frc[1] / self.m_total;
        self.v_com[2] = self.v_com[2] + dtb2 * self.frc[2] / self.m_total;

        // rotate the axis to next time
    }
    fn stage2(&mut self, dt: f32){

    }

}
