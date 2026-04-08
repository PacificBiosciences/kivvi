// Source code from trgt

#[derive(Debug, PartialEq, Clone)]
pub enum Color {
    Purple,
    Blue,
    Orange,
    Teal,
    Gray,
    LightGray,
    Black,
    Green,
    Pink,
    Yellow,
    Red,
    Khaki,
    PaleRed,
    PaleBlue,
    White,
    Grad(f64),
}

#[derive(Debug, PartialEq)]
pub enum Shape {
    Rect,
    HLine,
    VLine,
}

#[derive(Debug)]
pub struct PipeSeg {
    pub width: u32,
    pub color: Color,
    pub shape: Shape,
}

#[derive(Debug)]
pub struct Beta {
    pub value: f64,
    pub pos: usize,
}

#[derive(Debug)]
pub struct Pipe {
    pub segs: Vec<PipeSeg>,
    pub betas: Vec<Beta>,
    pub height: u32,
    pub outline: Vec<PipeSeg>,
    pub scale: Vec<(u32, Option<u32>)>,
}

pub struct Legend {
    pub labels: Vec<(String, Color)>,
    pub height: u32,
}

pub type PlotPanel = Vec<Pipe>;

pub struct PipePlot {
    pub panels: Vec<PlotPanel>,
    pub legend: Legend,
    // pub has_meth: bool,
}

pub fn encode_color(color: &Color) -> String {
    match color {
        Color::Purple => "#814ED1".to_string(),
        Color::Blue => "#1383C6".to_string(),
        Color::Orange => "#E16A2C".to_string(),
        Color::Teal => "#009CA2".to_string(),
        Color::Gray => "#BABABA".to_string(),
        Color::LightGray => "#D1D1D1".to_string(),
        Color::Black => "#000000".to_string(),
        Color::Pink => "#ED3981".to_string(),
        Color::Yellow => "#EFCD17".to_string(),
        Color::Green => "#009D4E".to_string(),
        Color::Red => "#E3371E".to_string(),
        Color::Khaki => "#F0E68C".to_string(),
        Color::PaleRed => "#FF4858".to_string(),
        Color::PaleBlue => "#46B2E8".to_string(),
        Color::White => "#FFFFFF".to_string(),
        Color::Grad(value) => get_gradient(*value),
    }
}

fn get_gradient(value: f64) -> String {
    let blue: (u8, u8, u8) = (0, 73, 255);
    let red: (u8, u8, u8) = (255, 0, 0);
    let mix_red = (blue.0 as f64 * (1.0 - value) + red.0 as f64 * value).round() as u8;
    let mix_green = (blue.1 as f64 * (1.0 - value) + red.1 as f64 * value).round() as u8;
    let mix_blue = (blue.2 as f64 * (1.0 - value) + red.2 as f64 * value).round() as u8;

    format!("#{:02X}{:02X}{:02X}", mix_red, mix_green, mix_blue)
}
