def multiqc_sgr_config():
    from multiqc import config

    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "scdynascope/stats": {"fn": "*scdynascope.*stats.json"},
        "scdynascope/umi_count": {
            "fn": "*scdynascope.umi_count.json",
        },
        "scdynascope/saturation": {
            "fn": "*scdynascope.saturation.json",
        },
        "scdynascope/median_gene": {
            "fn": "*scdynascope.median_gene.json",
        },
        "scdynascope/substitution": {"fn": "*scdynascope.substitution.json"},
        "scdynascope/tor": {"fn": "*scdynascope.tor.json"},
    }
    config.update_dict(config.sp, sgr_search_patterns)
