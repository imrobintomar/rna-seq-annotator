# RNA-Seq Data Annotator

A robust Python tool for annotating RNA sequencing data using standardized ontologies and classification systems. This tool provides automated annotation capabilities while maintaining strict adherence to established scientific ontologies.

## Features

- Support for multiple ontology formats (primarily OBO)
- Automated annotation using Gene Ontology (GO) and Sequence Ontology (SO)
- Comprehensive validation framework
- Detailed logging system
- Type-safe implementation
- CSV output support
- Extensible architecture for custom annotation rules
- Docker support for containerized deployment
- MongoDB integration for data persistence
- Redis caching for improved performance

## Prerequisites

### Local Installation
- Python 3.8+
- pandas
- typing
- logging

### Docker Deployment
- Docker 20.10+
- Docker Compose 2.0+
- At least 8GB RAM
- 20GB free disk space

## Installation

### Local Installation

1. Clone the repository:
```bash
git clone https://github.com/imrobintomar/rna-seq-annotator.git
cd rna-seq-annotator
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

### Docker Deployment

1. Clone the repository:
```bash
git clone https://github.com/imrobintomar/rna-seq-annotator.git
cd rna-seq-annotator
```

2. Create necessary directories:
```bash
mkdir -p data ontologies output logs
```

3. Start the services:
```bash
docker-compose up -d
```

4. Verify services are running:
```bash
docker-compose ps
```

## Usage

### Local Usage

```python
from rna_seq_annotator import RNASeqAnnotator

# Define ontology files
ontology_files = {
    'GO': 'path/to/gene_ontology.obo',
    'SO': 'path/to/sequence_ontology.obo'
}

# Initialize annotator
annotator = RNASeqAnnotator(ontology_files)

# Load and annotate sequences
sequence_data = pd.read_csv('your_rna_seq_data.csv')
annotated_data = annotator.annotate_sequence(
    sequence_data,
    required_ontologies=['GO', 'SO'],
    output_file='annotated_sequences.csv'
)
```

### Docker Usage

1. Place your input files:
```bash
cp your_rna_seq_data.csv data/
cp gene_ontology.obo ontologies/
cp sequence_ontology.obo ontologies/
```

2. Run annotation:
```bash
docker-compose exec rna-seq-annotator python annotate.py \
    --input /app/data/your_rna_seq_data.csv \
    --output /app/output/annotated_sequences.csv
```

3. Access results:
```bash
ls output/
```

## Docker Services

The Docker deployment includes three main services:

### RNA-Seq Annotator
- Main application container
- Resource limits: 2 CPUs, 4GB RAM
- Mounted volumes for data, ontologies, output, and logs

### MongoDB
- Persistent data storage
- Secure authentication enabled
- Port: 27017
- Resource limits: 1 CPU, 2GB RAM

### Redis
- Caching layer
- Password protected
- Port: 6379
- Resource limits: 0.5 CPU, 1GB RAM

## Environment Variables

### RNA-Seq Annotator
- `LOG_LEVEL`: Logging level (default: INFO)
- `MAX_WORKERS`: Number of worker processes (default: 4)
- `ONTOLOGY_DIR`: Directory for ontology files
- `OUTPUT_DIR`: Directory for output files

### MongoDB
- `MONGO_INITDB_ROOT_USERNAME`: MongoDB admin username
- `MONGO_INITDB_ROOT_PASSWORD`: MongoDB admin password
- `MONGO_INITDB_DATABASE`: Default database name

### Redis
- `REDIS_PASSWORD`: Redis authentication password

## Input Data Format

The tool expects RNA-seq data in CSV format with the following columns:
- sequence_id (required): Unique identifier for each sequence
- sequence (required): The RNA sequence data
- Additional metadata columns (optional)

Example:
```csv
sequence_id,sequence,tissue_type,condition
seq001,AUGCAUGCAUGC,liver,control
seq002,GCAUGCAUGCAU,liver,treated
```

## Output Format

The tool generates annotated data in CSV format with additional columns for each ontology:
- Original columns
- GO_annotation: Gene Ontology annotations
- SO_annotation: Sequence Ontology annotations
- validation_status: Validation results (if validation is performed)

## Customization

### Adding New Ontologies

1. Prepare your ontology file in OBO format
2. Add the ontology to the `ontology_files` dictionary:
```python
ontology_files = {
    'GO': 'gene_ontology.obo',
    'SO': 'sequence_ontology.obo',
    'YOUR_ONTOLOGY': 'your_ontology.obo'
}
```

### Docker Customization

Modify `docker-compose.yml` to adjust:
- Resource limits
- Volume mounts
- Environment variables
- Network configuration

## Logging

The tool provides comprehensive logging with different levels:
- INFO: General operation information
- WARNING: Non-critical issues
- ERROR: Critical issues that need attention
- DEBUG: Detailed information for debugging

Logs are available:
- Local deployment: `rna_seq_annotator.log`
- Docker deployment: `/app/logs/rna_seq_annotator.log`

## Troubleshooting

### Docker Issues
1. Check service status:
```bash
docker-compose ps
```

2. View logs:
```bash
docker-compose logs [service-name]
```

3. Restart services:
```bash
docker-compose restart [service-name]
```

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License 

## Citation

If you use this tool in your research, please cite:
```bibtex
@software{rna_seq_annotator,
  author = {Robin Tomar},
  title = {RNA-Seq Data Annotator},
  year = {2024},
  url = {https://github.com/imrobintomar/rna-seq-annotator}
}
```

## Support

For support:
- Open an issue in the GitHub repository
- Contact [itsrobintomar@gmail.com]
- Check the troubleshooting guide

## Acknowledgments

- Gene Ontology Consortium
- Sequence Ontology Project
- Contributors and maintainers of dependent packages
- Docker and container ecosystem contributors
