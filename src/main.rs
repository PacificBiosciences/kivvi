use kivvi::caller::{call_d4z4, call_kiv};
use kivvi::cli::{check_settings, get_raw_settings, Command, Settings};
use kivvi::util::{FULL_VERSION_PROGRAM, GIT_DESCRIBE};
use log::LevelFilter;

fn main() {
    // get cli settings
    let settings: Settings = get_raw_settings();
    let filter_level: LevelFilter = match settings.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    // setup logging
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        //.filter_level(LevelFilter::Warn)
        //.filter_module(env!("CARGO_PKG_NAME"), filter_level)
        .init();

    // run caller
    let cli_settings: Settings = check_settings(settings);

    let subcommand_name = match cli_settings.command {
        Command::Kiv2(_) => "KIV2",
        Command::D4z4(_) => "D4Z4",
    };

    log::info!(
        "Running {} ({}) [{}]",
        &*FULL_VERSION_PROGRAM,
        &*GIT_DESCRIBE,
        subcommand_name
    );
    match cli_settings.command {
        Command::Kiv2(_) => {
            if let Err(e) = call_kiv(cli_settings) {
                log::error!("{}", e);
                std::process::exit(1);
            }
        }
        Command::D4z4(_) => {
            if let Err(e) = call_d4z4(cli_settings) {
                log::error!("{}", e);
                std::process::exit(1);
            }
        }
    }
}
