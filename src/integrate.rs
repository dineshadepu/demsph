pub trait RK2 {
    fn initialize(&mut self, dt: f32);
    fn stage1(&mut self, dt: f32);
    fn stage2(&mut self, dt: f32);
}

pub fn integrate_initialize<T: RK2>(world: &mut Vec<&mut T>, dt: f32) {
    for entity in world {
        entity.initialize(dt)
    }
}

pub fn integrate_stage1<T: RK2>(world: &mut Vec<&mut T>, dt: f32) {
    for entity in world {
        entity.stage1(dt)
    }
}

pub fn integrate_stage2<T: RK2>(world: &mut Vec<&mut T>, dt: f32) {
    for entity in world {
        entity.stage2(dt)
    }
}
