use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};
use nalgebra as na;
use std::f64::consts::PI;

#[derive(Debug)]
struct Panel {
    start: na::Point2<f64>,
    end: na::Point2<f64>,
    control_point: na::Point2<f64>,
    normal: na::Vector2<f64>,
    length: f64,
}

impl Panel {
    fn new(start: na::Point2<f64>, end: na::Point2<f64>) -> Self {
        let midpoint = na::Point2::new(
            (start.x + end.x) * 0.5,
            (start.y + end.y) * 0.5,
        );

        let dx = end.x - start.x;
        let dy = end.y - start.y;
        let length = (dx * dx + dy * dy).sqrt();

        // Normal vector pointing outward
        let normal = na::Vector2::new(-dy / length, dx / length);

        Panel {
            start,
            end,
            control_point: midpoint,
            normal,
            length,
        }
    }
}

#[derive(Debug)]
struct Airfoil {
    panels: Vec<Panel>,
    num_panels: usize,
}

impl Airfoil {
    fn create_naca_4digit(code: &str, num_panels: usize) -> Self {
        let m = code[0..1].parse::<f64>().unwrap() / 100.0;
        let p = code[1..2].parse::<f64>().unwrap() / 10.0;
        let t = code[2..4].parse::<f64>().unwrap() / 100.0;

        let mut x_points = Vec::with_capacity(num_panels + 1);
        let theta_spacing = PI / num_panels as f64;

        // Generate cosine-spaced x coordinates
        for i in 0..=num_panels {
            let theta = i as f64 * theta_spacing;
            x_points.push(0.5 * (1.0 - f64::cos(theta)));
        }

        let mut upper_surface = Vec::new();
        let mut lower_surface = Vec::new();

        for &x in &x_points {
            // Calculate mean camber line
            let yc = if x < p {
                m * (x / p.powi(2)) * (2.0 * p - x)
            } else {
                m * ((1.0 - x) / (1.0 - p).powi(2)) * (1.0 + x - 2.0 * p)
            };

            // Calculate thickness distribution
            let yt = 5.0 * t * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2) +
                0.2843 * x.powi(3) - 0.1015 * x.powi(4));

            // Calculate surface slopes
            let dyc_dx = if x < p {
                2.0 * m / p.powi(2) * (p - x)
            } else {
                2.0 * m / (1.0 - p).powi(2) * (p - x)
            };

            let theta = f64::atan(dyc_dx);

            // Calculate upper and lower surface coordinates
            upper_surface.push(na::Point2::new(
                x - yt * f64::sin(theta),
                yc + yt * f64::cos(theta)
            ));

            lower_surface.push(na::Point2::new(
                x + yt * f64::sin(theta),
                yc - yt * f64::cos(theta)
            ));
        }

        // Combine points to form panels
        lower_surface.reverse();
        let surface_points = [upper_surface, lower_surface].concat();

        let panels = surface_points.windows(2)
            .map(|pts| Panel::new(pts[0], pts[1]))
            .collect();

        Airfoil {
            panels,
            num_panels,
        }
    }
}

#[derive(Debug)]
struct StreamlineData {
    points: Vec<[f64; 2]>,
    velocities: Vec<f64>,
}

impl StreamlineData {
    fn new() -> Self {
        Self {
            points: Vec::new(),
            velocities: Vec::new(),
        }
    }
}

#[derive(Debug)]
struct PressureDistribution {
    cp: Vec<f64>,              // Pressure coefficient at each panel
    x_positions: Vec<f64>,     // x/c positions along the surface
    is_upper_surface: Vec<bool>, // Whether point is on upper surface
}

#[derive(Debug)]
struct PerformanceMetrics {
    cl: f64,     // Lift coefficient
    cd: f64,     // Drag coefficient
    cm: f64,     // Moment coefficient (about quarter-chord)
}

#[derive(Debug)]
struct FlowField {
    freestream_velocity: f64,
    angle_of_attack: f64,
    density: f64,
    panels_solution: Vec<f64>,
}

impl FlowField {
    fn solve_panel_method(&mut self, airfoil: &Airfoil) {
        let n = airfoil.panels.len();
        let mut influence_matrix = na::DMatrix::zeros(n, n);
        let mut rhs = na::DVector::zeros(n);

        // Angle in radians
        let alpha = self.angle_of_attack * PI / 180.0;

        // Build influence coefficient matrix and RHS
        for i in 0..n {
            let control_point = &airfoil.panels[i].control_point;
            let normal = &airfoil.panels[i].normal;

            // RHS: freestream contribution
            rhs[i] = -(self.freestream_velocity *
                (f64::cos(alpha) * normal.x + f64::sin(alpha) * normal.y));

            // Influence coefficients
            for j in 0..n {
                let panel = &airfoil.panels[j];
                influence_matrix[(i, j)] = self.calculate_influence(
                    control_point,
                    panel,
                    normal
                );
            }
        }

        // Solve system using Gaussian elimination
        let decomposition = influence_matrix.lu();
        let solution = decomposition.solve(&rhs).unwrap();

        self.panels_solution = solution.as_slice().to_vec();
    }

    fn calculate_influence(&self, point: &na::Point2<f64>, panel: &Panel, normal: &na::Vector2<f64>) -> f64 {
        // Simple source panel method influence calculation
        let dx = point.x - panel.control_point.x;
        let dy = point.y - panel.control_point.y;
        let r_squared = dx * dx + dy * dy;

        if r_squared < 1e-10 {
            return 0.0;  // Avoid self-influence singularity
        }

        let r = r_squared.sqrt();
        let cos_theta = dx / r;
        let sin_theta = dy / r;

        (cos_theta * normal.x + sin_theta * normal.y) / (2.0 * PI * r)
    }

    fn calculate_velocity(&self, point: &na::Point2<f64>, airfoil: &Airfoil) -> na::Vector2<f64> {
        let alpha = self.angle_of_attack * PI / 180.0;

        // Freestream velocity components
        let mut vx = self.freestream_velocity * f64::cos(alpha);
        let mut vy = self.freestream_velocity * f64::sin(alpha);

        // Add influence of all panels
        for (i, panel) in airfoil.panels.iter().enumerate() {
            let strength = self.panels_solution[i];
            let dx = point.x - panel.control_point.x;
            let dy = point.y - panel.control_point.y;
            let r_squared = dx * dx + dy * dy;

            if r_squared < 1e-10 {
                continue;  // Skip if too close to avoid singularity
            }

            let r = r_squared.sqrt();
            let u = strength * dx / (2.0 * PI * r_squared);
            let v = strength * dy / (2.0 * PI * r_squared);

            vx += u;
            vy += v;
        }

        na::Vector2::new(vx, vy)
    }

    fn calculate_streamline(&self, start_point: na::Point2<f64>, airfoil: &Airfoil) -> StreamlineData {
        let mut streamline_data = StreamlineData::new();
        let mut current_point = start_point;
        streamline_data.points.push([current_point.x, current_point.y]);

        let dt = 0.01;  // Time step
        let steps = 200;  // Number of steps to calculate

        for _ in 0..steps {
            let velocity = self.calculate_velocity(&current_point, airfoil);
            let speed = (velocity.x * velocity.x + velocity.y * velocity.y).sqrt();

            streamline_data.velocities.push(speed);

            if speed < 1e-10 {
                break;  // Stop if velocity is too small
            }

            // Update position using RK4 integration
            let k1 = velocity;
            let p2 = current_point + (k1 * dt / 2.0);
            let k2 = self.calculate_velocity(&p2, airfoil);
            let p3 = current_point + (k2 * dt / 2.0);
            let k3 = self.calculate_velocity(&p3, airfoil);
            let p4 = current_point + (k3 * dt);
            let k4 = self.calculate_velocity(&p4, airfoil);

            let dp = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
            current_point += dp;

            streamline_data.points.push([current_point.x, current_point.y]);

            // Stop if streamline goes too far
            if current_point.x < -2.0 || current_point.x > 3.0 ||
                current_point.y < -2.0 || current_point.y > 2.0 {
                break;
            }
        }

        streamline_data
    }

    fn calculate_pressure_distribution(&self, airfoil: &Airfoil) -> PressureDistribution {
        let mut cp = Vec::with_capacity(airfoil.panels.len());
        let mut x_positions = Vec::with_capacity(airfoil.panels.len());
        let mut is_upper_surface = Vec::with_capacity(airfoil.panels.len());

        let v_inf = self.freestream_velocity;

        for (i, panel) in airfoil.panels.iter().enumerate() {
            // Calculate local velocity at panel control point
            let velocity = self.calculate_velocity(&panel.control_point, airfoil);
            let v_mag = (velocity.x * velocity.x + velocity.y * velocity.y).sqrt();

            // Calculate pressure coefficient
            let cp_val = 1.0 - (v_mag / v_inf).powi(2);

            cp.push(cp_val);
            x_positions.push(panel.control_point.x);
            // Determine if point is on upper surface (y > 0 for first half of panels)
            is_upper_surface.push(i < airfoil.panels.len() / 2);
        }

        PressureDistribution {
            cp,
            x_positions,
            is_upper_surface,
        }
    }

    fn calculate_performance_metrics(&self, airfoil: &Airfoil) -> PerformanceMetrics {
        let mut cl = 0.0;
        let mut cd = 0.0;
        let mut cm = 0.0;

        // Angle in radians
        let alpha = self.angle_of_attack * PI / 180.0;

        for (i, panel) in airfoil.panels.iter().enumerate() {
            let strength = self.panels_solution[i];
            let velocity = self.calculate_velocity(&panel.control_point, airfoil);
            let v_mag = (velocity.x * velocity.x + velocity.y * velocity.y).sqrt();

            // Calculate pressure coefficient
            let cp = 1.0 - (v_mag / self.freestream_velocity).powi(2);

            // Calculate force contribution from this panel
            let dx = panel.end.x - panel.start.x;
            let dy = panel.end.y - panel.start.y;
            let panel_length = (dx * dx + dy * dy).sqrt();

            // Force components in global coordinates
            let force_normal = -cp * panel_length;

            // Project forces to lift and drag directions
            let lift_contribution = force_normal *
                (panel.normal.x * f64::sin(alpha) - panel.normal.y * f64::cos(alpha));
            let drag_contribution = force_normal *
                (-panel.normal.x * f64::cos(alpha) - panel.normal.y * f64::sin(alpha));

            // Quarter-chord moment contribution
            let panel_moment_arm = panel.control_point.x - 0.25; // Relative to quarter-chord
            let moment_contribution = -cp * panel_length * panel_moment_arm;

            cl += lift_contribution;
            cd += drag_contribution;
            cm += moment_contribution;
        }

        // Normalize by freestream dynamic pressure and chord length
        let q_inf = 0.5 * self.density * self.freestream_velocity.powi(2);
        cl /= q_inf;
        cd /= q_inf;
        cm /= (q_inf);  // Additional chord length for moment coefficient

        PerformanceMetrics { cl, cd, cm }
    }
}

#[derive(Debug)]
struct VelocityFieldSettings {
    show_field: bool,
    grid_density: usize,    // Number of points in each direction
    scale_factor: f64,      // Scale factor for vector lengths
    min_velocity: f64,
    max_velocity: f64,
}

struct AirflowSimulator {
    airfoil: Airfoil,
    flow: FlowField,
    show_panels: bool,
    show_streamlines: bool,
    show_pressure: bool,
    naca_code: String,
    streamlines: Vec<StreamlineData>,
    min_velocity: f64,
    max_velocity: f64,
    pressure_distribution: Option<PressureDistribution>,
    performance_metrics: Option<PerformanceMetrics>,
    velocity_field_settings: VelocityFieldSettings,
    field_points: Option<Vec<na::Point2<f64>>>,
    field_velocities: Option<Vec<na::Vector2<f64>>>,
}

impl AirflowSimulator {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let airfoil = Airfoil::create_naca_4digit("0012", 50);

        Self {
            airfoil,
            flow: FlowField {
                freestream_velocity: 1.0,
                angle_of_attack: 5.0,
                density: 1.225,
                panels_solution: vec![],
            },
            show_panels: true,
            show_streamlines: true,
            show_pressure: false,
            naca_code: "0012".to_string(),
            streamlines: Vec::new(),
            min_velocity: f64::INFINITY,
            max_velocity: -f64::INFINITY,
            pressure_distribution: None,
            performance_metrics: None,
            velocity_field_settings: VelocityFieldSettings {
                show_field: false,
                grid_density: 20,     // 20x20 grid by default
                scale_factor: 0.1,    // scale factor to make vectors visible
                min_velocity: f64::INFINITY,
                max_velocity: -f64::INFINITY,
            },
            field_points: None,
            field_velocities: None,
        }
    }

    fn run_simulation(&mut self) {
        self.flow.solve_panel_method(&self.airfoil);
        self.calculate_streamlines();
        self.pressure_distribution = Some(self.flow.calculate_pressure_distribution(&self.airfoil));
        self.performance_metrics = Some(self.flow.calculate_performance_metrics(&self.airfoil));
        self.calculate_velocity_field();
    }

    fn calculate_streamlines(&mut self) {
        self.streamlines.clear();
        self.min_velocity = f64::INFINITY;
        self.max_velocity = -f64::INFINITY;

        // Generate streamline starting points
        let y_starts = (-10..=10).map(|i| i as f64 * 0.2);

        for y in y_starts {
            let start_point = na::Point2::new(-1.0, y);
            let streamline_data = self.flow.calculate_streamline(start_point, &self.airfoil);

            // Update min/max velocities
            if let Some(min) = streamline_data.velocities.iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {
                self.min_velocity = self.min_velocity.min(*min);
            }
            if let Some(max) = streamline_data.velocities.iter().max_by(|a, b| a.partial_cmp(b).unwrap()) {
                self.max_velocity = self.max_velocity.max(*max);
            }

            self.streamlines.push(streamline_data);
        }
    }

    fn get_airfoil_points(&self) -> PlotPoints {
        let points: Vec<[f64; 2]> = self.airfoil.panels.iter()
            .map(|panel| [panel.start.x, panel.start.y])
            .collect();

        PlotPoints::new(points)
    }

    fn velocity_to_color(&self, velocity: f64) -> egui::Color32 {
        // Normalize velocity to [0, 1] range
        let t = (velocity - self.min_velocity) / (self.max_velocity - self.min_velocity);

        // Create a color gradient from blue (slow) to red (fast)
        let r = (t * 255.0) as u8;
        let b = ((1.0 - t) * 255.0) as u8;
        let g = (((1.0 - t) * t * 4.0) * 255.0) as u8;  // Peak green in the middle

        egui::Color32::from_rgb(r, g, b)
    }

    fn plot_pressure_distribution(&self, plot_ui: &mut egui_plot::PlotUi) {
        if let Some(pressure) = &self.pressure_distribution {
            // Scale the pressure coefficients to be visible on the airfoil plot
            // We'll offset them from the airfoil surface by scaling and shifting
            let scale_factor = 0.2;

            // Plot upper surface pressure
            let upper_points: Vec<[f64; 2]> = pressure.x_positions.iter()
                .zip(pressure.cp.iter())
                .zip(pressure.is_upper_surface.iter())
                .filter(|((_, _), &is_upper)| is_upper)
                .map(|((x, cp), _)| {
                    // Place pressure distribution above the airfoil
                    [*x, cp * scale_factor + 0.1]
                })
                .collect();

            let upper_line = Line::new(upper_points)
                .color(egui::Color32::from_rgb(255, 100, 100))
                .width(2.0)
                .name("Upper Surface Cp");

            // Plot lower surface pressure
            let lower_points: Vec<[f64; 2]> = pressure.x_positions.iter()
                .zip(pressure.cp.iter())
                .zip(pressure.is_upper_surface.iter())
                .filter(|((_, _), &is_upper)| !is_upper)
                .map(|((x, cp), _)| {
                    // Place pressure distribution below the airfoil
                    [*x, cp * scale_factor - 0.1]
                })
                .collect();

            let lower_line = Line::new(lower_points)
                .color(egui::Color32::from_rgb(100, 100, 255))
                .width(2.0)
                .name("Lower Surface Cp");

            plot_ui.line(upper_line);
            plot_ui.line(lower_line);

            // Optionally add reference lines to show zero pressure
            let zero_upper = Line::new(vec![[0.0, 0.1], [1.0, 0.1]])
                .color(egui::Color32::from_rgba_premultiplied(100, 100, 100, 100))
                .width(1.0);
            let zero_lower = Line::new(vec![[0.0, -0.1], [1.0, -0.1]])
                .color(egui::Color32::from_rgba_premultiplied(100, 100, 100, 100))
                .width(1.0);

            plot_ui.line(zero_upper);
            plot_ui.line(zero_lower);
        }
    }

    fn show_performance_metrics(&self, ui: &mut egui::Ui) {
        if let Some(metrics) = &self.performance_metrics {
            ui.heading("Performance Metrics");

            // Create a frame with a light background
            egui::Frame::none()
                .fill(ui.visuals().extreme_bg_color)
                .show(ui, |ui| {
                    ui.vertical(|ui| {
                        // Display metrics with consistent formatting
                        ui.label(format!("Lift Coefficient (CL): {:.3}", metrics.cl));
                        ui.label(format!("Drag Coefficient (CD): {:.4}", metrics.cd));
                        ui.label(format!("Moment Coefficient (CM): {:.3}", metrics.cm));

                        // Calculate and display lift-to-drag ratio
                        let l_d_ratio = metrics.cl / metrics.cd;
                        ui.label(format!("Lift-to-Drag Ratio (L/D): {:.1}", l_d_ratio));
                    });
                });
        }
    }

    fn calculate_velocity_field(&mut self) {
        let settings = &mut self.velocity_field_settings;
        settings.min_velocity = f64::INFINITY;
        settings.max_velocity = -f64::INFINITY;

        // Create a grid of points
        let x_min = -0.5;
        let x_max = 1.5;
        let y_min = -1.0;
        let y_max = 1.0;

        let dx = (x_max - x_min) / settings.grid_density as f64;
        let dy = (y_max - y_min) / settings.grid_density as f64;

        let mut field_points = Vec::new();
        let mut field_velocities = Vec::new();
        
        for i in 0..=settings.grid_density {
            for j in 0..=settings.grid_density {
                let x = x_min + dx * i as f64;
                let y = y_min + dy * j as f64;
                let point = na::Point2::new(x, y);
                
                let mut is_near_point: bool = false;

                for panel in &self.airfoil.panels {
                    let dx = point.x - panel.control_point.x;
                    let dy = point.y - panel.control_point.y;
                    let distance_squared = dx * dx + dy * dy;

                    if distance_squared < 0.01 {
                        is_near_point = true;
                    }
                }

                // Skip points inside or very close to the airfoil
                if is_near_point {
                    continue;
                }

                let velocity = self.flow.calculate_velocity(&point, &self.airfoil);
                let speed = (velocity.x * velocity.x + velocity.y * velocity.y).sqrt();

                settings.min_velocity = settings.min_velocity.min(speed);
                settings.max_velocity = settings.max_velocity.max(speed);

                field_points.push(point);
                field_velocities.push(velocity);
            }
        }

        self.field_points = Some(field_points);
        self.field_velocities = Some(field_velocities);
    }

    fn is_point_near_airfoil(&self, point: &na::Point2<f64>) -> bool {
        for panel in &self.airfoil.panels {
            let dx = point.x - panel.control_point.x;
            let dy = point.y - panel.control_point.y;
            let distance_squared = dx * dx + dy * dy;

            if distance_squared < 0.01 { // threshold
                return true;
            }
        }
        false
    }

    fn draw_velocity_field(&self, plot_ui: &mut egui_plot::PlotUi) {
        if let (Some(points), Some(velocities)) = (&self.field_points, &self.field_velocities) {
            let settings = &self.velocity_field_settings;

            for (point, velocity) in points.iter().zip(velocities.iter()) {
                let speed = (velocity.x * velocity.x + velocity.y * velocity.y).sqrt();

                // Normalize velocity to [0, 1] range for color
                let t = (speed - settings.min_velocity) /
                    (settings.max_velocity - settings.min_velocity);

                // Create color gradient from blue (slow) to red (fast)
                let color = self.velocity_to_color(speed);

                // Calculate arrow end point
                let scaled_dx = velocity.x * settings.scale_factor / settings.max_velocity;
                let scaled_dy = velocity.y * settings.scale_factor / settings.max_velocity;

                // Draw arrow line
                let arrow = Line::new(vec![
                    [point.x, point.y],
                    [point.x + scaled_dx, point.y + scaled_dy],
                ])
                    .color(color)
                    .width(1.0);

                plot_ui.line(arrow);

                // Draw arrow head
                // TODO: find a better less ghetto to do this
                let head_size = 0.02;
                let angle = f64::atan2(scaled_dy, scaled_dx);
                let head_angle = 0.5;  // Radians

                let arrow_head = Line::new(vec![
                    [point.x + scaled_dx, point.y + scaled_dy],
                    [
                        point.x + scaled_dx - head_size * f64::cos(angle + head_angle),
                        point.y + scaled_dy - head_size * f64::sin(angle + head_angle),
                    ],
                    [
                        point.x + scaled_dx - head_size * f64::cos(angle - head_angle),
                        point.y + scaled_dy - head_size * f64::sin(angle - head_angle),
                    ],
                    [point.x + scaled_dx, point.y + scaled_dy],
                ])
                    .color(color)
                    .width(1.0);

                plot_ui.line(arrow_head);
            }
        }
    }
}

impl eframe::App for AirflowSimulator {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        
        egui::SidePanel::left("controls").show(ctx, |ui| {
            ui.heading("Simulation Controls");

            ui.horizontal(|ui| {
                ui.label("NACA:");
                if ui.text_edit_singleline(&mut self.naca_code).changed() {
                    if self.naca_code.len() == 4 {
                        self.airfoil = Airfoil::create_naca_4digit(&self.naca_code, 50);
                        self.streamlines.clear();
                    }
                }
            });

            ui.add(egui::Slider::new(&mut self.flow.angle_of_attack, -10.0..=10.0)
                .text("Angle of Attack (°)"));

            ui.add(egui::Slider::new(&mut self.flow.freestream_velocity, 0.1..=200.0)
                .text("Velocity (m/s)"));

            ui.checkbox(&mut self.show_panels, "Show Panels");
            ui.checkbox(&mut self.show_streamlines, "Show Streamlines");
            ui.checkbox(&mut self.show_pressure, "Show Pressure");

            if ui.button("Run Simulation").clicked() {
                self.run_simulation();
            }

            ui.separator();
            if let Some(_) = &self.performance_metrics {
                self.show_performance_metrics(ui);
            } else {
                ui.label("Run simulation to see performance metrics");
            }

            ui.separator();
            ui.heading("Velocity Field");
            ui.checkbox(&mut self.velocity_field_settings.show_field, "Show Velocity Vectors");

            if self.velocity_field_settings.show_field {
                ui.add(egui::Slider::new(&mut self.velocity_field_settings.grid_density, 10..=40)
                    .text("Grid Density"));
                ui.add(egui::Slider::new(&mut self.velocity_field_settings.scale_factor, 0.05..=0.3)
                    .text("Vector Scale"));
            }
            
            ui.vertical(|ui| {
                ui.label(format!("Velocity Scale (m/s):"));
                ui.label(format!("Max: {:.2}", self.max_velocity));
                ui.label(format!("Min: {:.2}", self.min_velocity));
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Airflow Visualization");
            
            let plot = Plot::new("airfoil_plot")
                .view_aspect(1.0)
                .include_x(-1.5)
                .include_x(2.0)
                .include_y(-1.5)
                .include_y(1.5);

            plot.show(ui, |plot_ui| {
                // Draw airfoil surface
                let airfoil_line = Line::new(self.get_airfoil_points())
                    .color(egui::Color32::DARK_RED)
                    .width(2.0);
                plot_ui.line(airfoil_line);

                // Draw streamlines
                if self.show_streamlines && !self.streamlines.is_empty() {
                    for streamline in &self.streamlines {
                        // Draw line segments with colors based on velocity
                        for i in 1..streamline.points.len() {
                            let velocity = streamline.velocities[i - 1];
                            let color = self.velocity_to_color(velocity);

                            let segment = Line::new(vec![
                                streamline.points[i - 1],
                                streamline.points[i],
                            ])
                                .color(color)
                                .width(1.5);

                            plot_ui.line(segment);
                        }
                    }
                }

                // Draw panel normal vectors
                if self.show_panels && !self.flow.panels_solution.is_empty() {
                    for (i, panel) in self.airfoil.panels.iter().enumerate() {
                        let strength = self.flow.panels_solution[i];
                        let normal_line = Line::new(vec![
                            [panel.control_point.x, panel.control_point.y],
                            [
                                panel.control_point.x + panel.normal.x * strength * 0.1,
                                panel.control_point.y + panel.normal.y * strength * 0.1,
                            ]
                        ])
                            .color(egui::Color32::RED)
                            .width(1.0);
                        plot_ui.line(normal_line);
                    }
                }

                if self.show_pressure {
                    self.plot_pressure_distribution(plot_ui);
                }

                if self.velocity_field_settings.show_field {
                    self.draw_velocity_field(plot_ui);
                }
            });

            if self.show_pressure {
                ui.horizontal(|ui| {
                    ui.label("Pressure Distribution:");
                    ui.label(egui::RichText::new("— Upper Surface")
                        .color(egui::Color32::from_rgb(255, 100, 100)));
                    ui.label(egui::RichText::new("— Lower Surface")
                        .color(egui::Color32::from_rgb(100, 100, 255)));
                });
            }

            if self.velocity_field_settings.show_field {
                ui.horizontal(|ui| {
                    ui.label("Velocity Magnitude:");
                    ui.colored_label(
                        egui::Color32::from_rgb(0, 0, 255),
                        format!("{:.1} m/s", self.velocity_field_settings.min_velocity)
                    );
                    ui.label("→");
                    ui.colored_label(
                        egui::Color32::from_rgb(255, 0, 0),
                        format!("{:.1} m/s", self.velocity_field_settings.max_velocity)
                    );
                });
            }
        });
    }
}

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Airflow Simulator",
        options,
        Box::new(|cc| Box::new(AirflowSimulator::new(cc)))
    )
}