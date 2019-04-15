use libc::{c_float, size_t, c_double};
use splines::{Interpolation, Key, Spline};
use bspline;
use std::mem;
use std::slice;
use std::vec::Vec;
use peroxide::numerical::spline;
use peroxide::Polynomial;
pub struct Bspline {
    spline: Spline<f32>,
}

#[no_mangle]
pub extern "C" fn new_spline(
    x: *const c_float,
    y: *const c_float,
    len: size_t,
    itype: u8,
) -> *mut Bspline {
    let mut keys: Vec<Key<f32>> = Vec::new();
    let newx = unsafe { slice::from_raw_parts(x, len as usize) };
    let newy = unsafe { slice::from_raw_parts(y, len as usize) };
    let intp = match itype {
        1 => Interpolation::Linear,
        2 => Interpolation::Cosine,
        3 => Interpolation::CatmullRom,
        _ => Interpolation::default(),
    };
    for i in 0..len as usize {
        keys.push(Key::new(newx[i], newy[i], intp));
    }

    let mut spline = Spline::from_vec(keys);
    //    println!("{:?}", spline);
    Box::into_raw(Box::new(Bspline { spline }))
}

#[no_mangle]
pub extern "C" fn sample_spline(spline: *mut Bspline, t: c_float) -> c_float {
    let spl = unsafe { &mut *spline };
    match spl.spline.sample(t) {
        Some(f) => f,
        _ => 0.0,
    }
}

#[no_mangle]
pub extern "C" fn free_spline(ptr: *mut Bspline) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        Box::from_raw(ptr);
    }
}


pub struct BSpline {
    spline: bspline::BSpline<f32>,
}

#[repr(C)]
pub struct Tuple {
    x:c_float,
    y:c_float,
}
// Conversion functions
impl From<(f32, f32)> for Tuple {
    fn from(tup: (f32, f32)) -> Tuple {
        Tuple { x: tup.0, y: tup.1 }
    }
}

impl From<Tuple> for (f32, f32) {
    fn from(tup: Tuple) -> (f32, f32) {
        (tup.x, tup.y)
    }
}

#[no_mangle]
pub extern "C" fn new_bspline(points:*const c_float,knots:*const c_float,lpts:size_t,lkts:size_t) -> *mut BSpline
{
    let degree = 3;
    let _pts = unsafe { slice::from_raw_parts(points, lpts as usize)};
    let _kts = unsafe { slice::from_raw_parts(knots, lkts as usize)};
    let vpoints = Vec::from(_pts);
    let vknots = Vec::from(_kts);
    let spline = bspline::BSpline::new(degree,vpoints,vknots);
    Box::into_raw(Box::new(BSpline{spline}))
}

#[no_mangle]
pub extern fn knot_domain(spline:*mut BSpline) -> Tuple
{
    let _spline = unsafe{ &mut *spline};
    _spline.spline.knot_domain().into()
}
#[no_mangle]
pub extern fn sample(spline:*mut BSpline,t:c_float) -> c_float
{
    let _spline = unsafe{ &mut *spline};
    _spline.spline.point(t)
}

pub struct pBSpline
{
    spline:Vec<Polynomial>,
}

#[no_mangle]
pub extern "C" fn new_bspline_peroxide(x:*const c_double,y:*const c_double,lpts:size_t,lkts:size_t) -> *mut pBSpline
{
    let _x = unsafe { slice::from_raw_parts(x, lpts as usize)};
    let _y = unsafe { slice::from_raw_parts(y, lkts as usize)};
    let vx = Vec::from(_x);
    let vy = Vec::from(_y);
    let spline = spline::cubic_spline(vx, vy);
    let npolynomials = spline.len() as usize;
    Box::into_raw(Box::new(pBSpline{spline}))
}

#[no_mangle]
pub extern fn eval_poly(spline:*mut pBSpline,outp:*mut f64,lout:size_t,t:c_double)
{
    let _spline = unsafe{ &mut *spline};
    let pvec = &_spline.spline;

    let mut retv:[f64;100] = [0.0;100];

    for (i,p) in _spline.spline.iter().enumerate()
    {

        retv[i]=p.eval(t);
        println!("{}",retv[i]);
    }
    unsafe {
    std::slice::from_raw_parts_mut(outp, lout)
        .copy_from_slice(&retv);
    }
}
// #[no_mangle]
// pub extern fn sample(spline:*mut BSpline,t:c_float) -> c_float
// {
//     let _spline = unsafe{ &mut *spline};
//     _spline.spline.point(t)
// }
