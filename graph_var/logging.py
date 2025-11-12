"""
Custom logging functionality for pantree with memory tracking.
"""
import logging
import psutil
from datetime import datetime
from typing import Optional


class MemoryFormatter(logging.Formatter):
    """Custom formatter that adds memory usage to log messages."""
    
    def format(self, record):
        # Measure memory in MB
        memory_mb = psutil.Process().memory_info().rss / 1024 / 1024
        
        # Add memory info to the record
        record.memory_mb = f"{memory_mb:.2f}"
        
        return super().format(record)


def setup_logger(name: str = "pantree", 
                 log_path: Optional[str] = None,
                 verbose: bool = False) -> logging.Logger:
    """Set up a logger with optional file output and console output.
    
    Args:
        name: Logger name
        log_path: Optional path to log file
        verbose: If True, also log to console
        
    Returns:
        Configured Logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers to avoid duplicates
    logger.handlers.clear()
    
    # Create formatter with timestamp and memory
    formatter = MemoryFormatter('[%(asctime)s] [%(memory_mb)s MB] %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S')
    
    # Add file handler if log_path provided
    if log_path:
        file_handler = logging.FileHandler(log_path, mode='a')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    # Add console handler if verbose
    if verbose:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    return logger
