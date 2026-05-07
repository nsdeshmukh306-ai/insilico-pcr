# Dashboard Screenshots

This directory holds screenshots for the README. They must be taken from the live dashboard.

## How to Generate

```bash
python insilico_pcr/webapp/run.py
# → http://localhost:8765
# Click "Load demo", then take screenshots
```

## Naming Convention & Recommended Dimensions

| Filename | Panel to capture | Width × Height |
|---|---|---|
| `dashboard_overview.png` | Full page — run panel + stats row visible | 1440 × 900 |
| `alignment_viewer.png` | Amplicons & Alignment tab — click row 1 to open alignment | 1440 × 800 |
| `thermodynamics_panel.png` | Thermodynamics tab — Tm chart + radar visible | 1440 × 800 |
| `offtarget_explorer.png` | Off-target Explorer tab | 1440 × 700 |
| `genome_track.png` | Genome Overview tab — track + histogram | 1440 × 700 |
| `live_parameters.png` | Live Parameters tab — sliders visible | 1440 × 700 |
| `primer_quality.png` | Primer Quality tab — both primer cards | 1440 × 700 |

## Screenshot Tips

- Use browser zoom at 90% for more content per screenshot
- Load the demo data (not custom) for consistent, reproducible screenshots
- Use Chrome or Firefox DevTools device emulation at 1440px width for consistent sizing
- Crop to the content area only (exclude browser chrome)

## Animated GIF

An animated GIF (`../gifs/demo_walkthrough.gif`) showing the demo load + tab navigation is planned.
Create with [LICEcap](https://www.cockos.com/licecap/) or OBS Studio.
Recommended size: 1200 × 700 px, 30 fps, 15–20 seconds, < 10 MB.
