#!/usr/bin/env bash
# Recorta o corredor Benfica-Centro-Aldeota de um extrato .osm.pbf do Ceará
# e mantém só as vias por onde passam carros (mesmo filtro "drive" do osmnx).
#
# Uso (precisa do osmium-tool: sudo apt install osmium-tool):
#   curl -O https://download.openstreetmap.fr/extracts/south-america/brazil/northeast/ceara-latest.osm.pbf
#   bash codigo/recortar_osm.sh ceara-latest.osm.pbf dados/osm/corredor.osm
#   python codigo/otimizacao_antes_depois.py --osm-xml dados/osm/corredor.osm
set -euo pipefail
ENTRADA="$1"
SAIDA="$2"
TMP=$(mktemp -d)

# 1. Retângulo (oeste,sul,leste,norte): o mesmo BBOX do script Python
osmium extract -b -38.548,-3.752,-38.490,-3.718 --strategy complete_ways \
    "$ENTRADA" -o "$TMP/recorte.osm.pbf"

# 2. Só vias de automóvel (sem calçadas, ciclovias, vias de serviço...)
osmium tags-filter "$TMP/recorte.osm.pbf" \
    w/highway=motorway,motorway_link,trunk,trunk_link,primary,primary_link,secondary,secondary_link,tertiary,tertiary_link,unclassified,residential,living_street,road,busway \
    -o "$TMP/vias.osm.pbf"

# 3. Tira vias privadas ou proibidas para carros
osmium tags-filter -i "$TMP/vias.osm.pbf" \
    w/access=private,no w/motor_vehicle=no w/motorcar=no \
    -o "$SAIDA" --overwrite

rm -r "$TMP"
echo "Recorte salvo em $SAIDA"
